#define _USE_MATH_DEFINES
#include <GL/freeglut.h>
#include <vector>
#include <glm/glm.hpp>

const int WIDTH=800;
const int HEIGHT=600;
const int HALF_WIDTH = WIDTH>>1;
const int HALF_HEIGHT = HEIGHT>>1;

struct DistanceConstraint {	int p1, p2; float rest_length, k_prime; };
struct BendingConstraint {	int p1, p2, p3;	float rest_length,  w, k_prime;};

std::vector<GLushort> indices;
std::vector<DistanceConstraint> d_constraints;

std::vector<BendingConstraint> b_constraints;
//particle system
std::vector<glm::vec3> X; //position
std::vector<glm::vec3> tmp_X; //predicted position
std::vector<glm::vec3> V; //velocity
std::vector<glm::vec3> F;
std::vector<float> W; //inverse particle mass 
std::vector<glm::vec3> Ri; //Ri = Xi-Xcm  
int mouseX=0, mouseY=0; 
int state =1 ;  

const size_t solver_iterations = 10; //number of solver iterations per step. PBD  

const float global_damping=0.9985f;  //global simulation damping
const float kBend = 0.5f; 
const float kStretch = 0.5f; 
const float kDamp = 0.001f;
const glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);  

const float mass = 1.f;
const float dt = 1.0f/60;
const int N = 10;
const float EPSILON = 0.0000001f;

int intersected = -1;

void AddDistanceConstraint(int a, int b, float k) {
	DistanceConstraint c;
	c.p1=a;
	c.p2=b; 
	c.k_prime = 1.0f-pow((1.0f-k), 1.0f/solver_iterations);   
	
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	 
	glm::vec3 deltaP = X[c.p1]-X[c.p2];
	c.rest_length = glm::length(deltaP);

	d_constraints.push_back(c);
}

void AddBendingConstraint(int pa, int pb, int pc, float k) {
	BendingConstraint c;
	c.p1=pa;
	c.p2=pb; 
	c.p3=pc;
	c.w = W[pa] + W[pb] + 2*W[pc]; 
	glm::vec3 center = 0.3333f * (X[pa] + X[pb] + X[pc]);
 	c.rest_length = glm::length(X[pc]-center);
	c.k_prime = 1.0f-pow((1.0f-k), 1.0f/solver_iterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints.push_back(c);
}

float map(float value, float istart, float istop, float ostart, float ostop) {
	return ostart + (ostop - ostart) * ((value - istart) / (istop - istart));
}

void OnInit() {

	glClearColor(200.0f/255, 200.0f/255, 200.0f/255,1);
	glPointSize(6);
	glLineWidth(5);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	X.resize(N);
	V.resize(N);
	F.resize(N);
	W.resize(N); 
	tmp_X.resize(N);
	Ri.resize(N);
   
	memset(&F[0],0,sizeof(glm::vec3)*N);
	memset(&V[0],0,sizeof(glm::vec3)*N);
	memset(&Ri[0],0,sizeof(glm::vec3)*N);
	memset(&tmp_X[0],0,sizeof(glm::vec3)*N);

    float dL = HEIGHT-float(HEIGHT)/N;
	for(int i=0; i<N; ++i) {    
		X[i] = glm::vec3( map(float(i),0.0f,float(N),float(HALF_WIDTH),float(WIDTH)), dL,0.0f); 
		W[i] = (i==0)?0.0f:(1.0f/mass); 
	} 
  
	for(int i=0; i<(N-1); ++i) {    
		AddDistanceConstraint(i, i+1, kStretch);
	}
  
	for(int i=0; i<(N-2); ++i) {    
		AddBendingConstraint(i, i+1, i+2, kBend);
	} 
}

void ComputeForces( ) {
	size_t i=0;
	
	for(i=0;i<N;i++) {
		F[i] = glm::vec3(0.0f); 
		 
		//add gravity force
		if(W[i]!=0)		 
			F[i] += gravity; 
	}	
} 
 
void IntegrateExplicitWithDamping(float deltaTime) { 
	size_t i=0;
 
	glm::vec3 Xcm = glm::vec3(0.0f);
	glm::vec3 Vcm = glm::vec3(0.0f);
	float sumM = 0;
	for(i=0;i<N;i++) {

		V[i] *= global_damping; //global velocity dampening !!!		
		V[i] = V[i] + (F[i]*deltaTime)*W[i]; 	 					
		
		//calculate the center of mass's position 
		//and velocity for damping calc
		Xcm += (X[i]*mass);
		Vcm += (V[i]*mass);
		sumM += mass;
	}
	Xcm /= sumM; 
	Vcm /= sumM; 

	glm::mat3 I = glm::mat3(1.0f);
	glm::vec3 L = glm::vec3(0.0f);
	glm::vec3 w = glm::vec3(0.0f);//angular velocity
	
	
	for(i=0;i<N;i++) { 
		Ri[i] = (X[i] - Xcm);	
		
		L += glm::cross(Ri[i],mass*V[i]); 

		//thanks to DevO for pointing this and these notes really helped.
		//http://www.sccg.sk/~onderik/phd/ca2010/ca10_lesson11.pdf

		glm::mat3 tmp = glm::mat3(0,-Ri[i].z,  Ri[i].y, 
							 Ri[i].z,       0,-Ri[i].x,
							 -Ri[i].y,Ri[i].x,        0);
		I +=(tmp*glm::transpose(tmp))*mass;
	} 
	
	w = glm::inverse(I)*L;
	
	//apply center of mass damping
	for(i=0;i<N;i++) {
		glm::vec3 delVi = Vcm + (glm::cross(w,Ri[i])-V[i]);		
		V[i] += kDamp*delVi;
	}

	//calculate predicted position
	for(i=0;i<N;i++) {
		if(W[i] <= 0.0) { 
			tmp_X[i] = X[i]; //fixed points
		} else {
			tmp_X[i] = X[i] + (V[i]*deltaTime);				 
		}
	} 
}
 
void Integrate(float deltaTime) {	
	float inv_dt = 1.0f/deltaTime;
	size_t i=0; 

	for(i=0;i<N;i++) {	
		V[i] = (tmp_X[i] - X[i])*inv_dt;		
		X[i] = tmp_X[i];		 
	}
}

void UpdateDistanceConstraint(int i) {

	DistanceConstraint c = d_constraints[i];
	glm::vec3 dir = tmp_X[c.p1] - tmp_X[c.p2];
	float len = glm::length(dir); 
	if(len <= EPSILON) 
		return;
	
	float w1 = W[c.p1];
	float w2 = W[c.p2];
	float invMass = w1+ w2; 
	if(invMass <= EPSILON) 
		return;
	 
	glm::vec3 dP = (1.0f/invMass) * (len-c.rest_length ) * dir/len * c.k_prime;
	if(w1 > 0.0)
		tmp_X[c.p1] -= dP*w1;

	if(w2 > 0.0)
		tmp_X[c.p2] += dP*w2;	
}

void UpdateBendingConstraint(int index) { 
	BendingConstraint c = b_constraints[index]; 

 	//Using the paper suggested by DevO
	//http://image.diku.dk/kenny/download/kelager.niebe.ea10.pdf
 
	//global_k is a percentage of the global dampening constant 
	float global_k = global_damping *0.01f; 
	glm::vec3 center = 0.3333f * (tmp_X[c.p1] + tmp_X[c.p2] + tmp_X[c.p3]);
	glm::vec3 dir_center = tmp_X[c.p3]-center;
	float dist_center = glm::length(dir_center);

	float diff = 1.0f - ((global_k + c.rest_length) / dist_center);
	glm::vec3 dir_force = dir_center * diff;
	glm::vec3 fa = c.k_prime * ((2.0f*W[c.p1])/c.w) * dir_force;
	glm::vec3 fb = c.k_prime * ((2.0f*W[c.p2])/c.w) * dir_force;
	glm::vec3 fc = -c.k_prime * ((4.0f*W[c.p3])/c.w) * dir_force;

	if(W[c.p1] > 0.0)  {
		tmp_X[c.p1] += fa;
	}
	if(W[c.p2] > 0.0) {
		tmp_X[c.p2] += fb;
	}
	if(W[c.p3] > 0.0) {
		tmp_X[c.p3] += fc;
	} 
}
 
void GroundCollision() 
{
	for(size_t i=0;i<N;i++) {	
		tmp_X[i].x = std::min<float>(tmp_X[i].x+2.0f, float(WIDTH));
		tmp_X[i].x = std::max<float>(tmp_X[i].x-2.0f, 0.0f);
     
		tmp_X[i].y = std::min<float>(tmp_X[i].y+2.0f, float(HEIGHT));
		tmp_X[i].y = std::max<float>(tmp_X[i].y-2.0f, 0.0f);
	}
}

void UpdateExternalConstraints() {
	GroundCollision();
}
//----------------------------------------------------------------------------------------------------
void UpdateInternalConstraints(float deltaTime) {
	size_t i=0;
 
	//printf(" UpdateInternalConstraints \n ");
	for (size_t si=0;si<solver_iterations;++si) {
		for(i=0;i<d_constraints.size();i++) {
			UpdateDistanceConstraint(i);
		} 
		for(i=0;i<b_constraints.size();i++) {
			UpdateBendingConstraint(i);
		} 
	}
}
void StepPhysics() {
	ComputeForces();
	IntegrateExplicitWithDamping(dt);
	UpdateInternalConstraints(dt);
	UpdateExternalConstraints();
	Integrate(dt);
}

void OnRender() {
	glClear(GL_COLOR_BUFFER_BIT); 
	StepPhysics();

	//lines using OpenGL
	glBegin(GL_LINES);
	glColor4f(0.5f, 0.5f, 0.5f, 120.0f/255);  
	for(int i=0; i<N-1; ++i) {
		glVertex2f(X[i].x, X[i].y);
		glVertex2f(X[i+1].x, X[i+1].y);
	}
	glEnd();
  
	glBegin(GL_POINTS);
	for(int i=0; i<N; ++i) {
		if(intersected==i)
			glColor3f(0,1,0);
		else
			glColor3f(0,0,0);
		glVertex2f(X[i].x, X[i].y);
	}
	glEnd();
	glutSwapBuffers();
}

void OnResize(int nw, int nh) {
	glViewport(0,0,nw,nh);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0,nw, 0, nh);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void OnIdle() {
	glutPostRedisplay();
}

void OnMouseDown(int button, int s, int x, int y)
{
	if (s == GLUT_DOWN) 
	{
		intersected = -1;
		mouseX = x;
		mouseY = HEIGHT-y;
		size_t i=0;
		glm::vec3 pt(x,y,0);
		for(i=0;i<N;i++) {			 
			float distX = mouseX-X[i].x;
			float distY = mouseY-X[i].y;
			
			if( (distX*distX + distY*distY) < 400) {
			    intersected = i;
				printf("Intersected at %d\n",i);
				break;
			}
		}
	}	

	if(button == GLUT_MIDDLE_BUTTON)
		state = 0;
	else
		state = 1;

	if(s==GLUT_UP) {
		intersected = -1;
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}

void OnMouseMove(int x, int y)
{
	if(intersected != -1) {
		X[intersected].x = float(mouseX);
		X[intersected].y = float(mouseY);
		V[intersected].x = 0;
		V[intersected].y = 0;	
	}
	
	mouseX = x; 
 	mouseY = HEIGHT-y; 
	
	glutPostRedisplay(); 
}



int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("PBD in 2D");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnResize);
	glutIdleFunc(OnIdle);
	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);

	OnInit();

	glutMainLoop();
	
	return 0;
}
