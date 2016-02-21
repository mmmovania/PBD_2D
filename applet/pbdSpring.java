import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class pbdSpring extends PApplet {

//A simple 2D implementation of the position based dynamics paper 
//by Dr. Muhammad Mobeen Movania. This demo shows how to model
//a 2D chain.
//
//You may click any point and drag to move it around

class DistanceConstraint {     
  public int p1, p2;     
  public float rest_length; 
  public float k_prime; 
};

class BendingConstraint {      
  public int p1, p2; 
  public float rest_length,  w,  k_prime;
};

final int solver_iterations= 10;//Total Gauss Seidel iterations 

PVector  X[];     //particle position
PVector tmp_X[];  //positional correction (p) in PBD paper
PVector  F[];     //Net force
PVector  V[];     //particle velocity
float    W[];     //particle inverse masses
float   phi0[];   //dihedral angle 
PVector Ri[];     //Ri=Xi-Xcm
PVector gravity=new PVector(0,9.8f);

ArrayList<DistanceConstraint> d_constraints;
ArrayList<BendingConstraint> b_constraints;

final float dt=1.0f/10.0f;  //time step (delT)
final int   N=10;         //total number of particles
final float mass=1;       //mass of each particle
final float global_damping=0.995f;  //global simulation damping
final float kBend = 0.5f;     //bending constraint coefficient
final float kStretch = 0.5f;  //distance constraint coefficient 
final float kDamp=0.001f;     //velocity damping coefficient
final float rest_length=50;  //constraint rest length
int intersected = -1;        //intersection flag

public void AddDistanceConstraint(int a, int b, float k) {
  DistanceConstraint c=new DistanceConstraint();
  c.p1=a;
  c.p2=b; 
  c.k_prime = 1.0f-pow((1.0f-k), 1.0f/solver_iterations);

  if(c.k_prime>1.0f)  
    c.k_prime = 1.0f;
         
  PVector deltaP = PVector.sub(X[c.p1],X[c.p2]);
  c.rest_length = deltaP.mag();
  d_constraints.add(c);
} 

public void AddBendingConstraint(int pa, int pb, float k) {
  BendingConstraint c=new BendingConstraint();
  c.p1=pa;
  c.p2=pb; 
        
  c.w = W[pa] + W[pb]; 
  PVector center = PVector.mult( PVector.add(X[pa],X[pb]), 0.5f);
  c.rest_length = PVector.sub(X[pa],center).mag();
  
  c.k_prime = 1.0f-pow((1.0f-k), 1.0f/solver_iterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);
  if(c.k_prime>1.0f) 
    c.k_prime = 1.0f;
  b_constraints.add(c);
}
 

public void setup(){
  size(800,600);
  
  X=new PVector[N];
  V=new PVector[N];
  F=new PVector[N];
  W=new float[N];
  phi0=new float[N];
  tmp_X=new PVector[N];
  Ri = new PVector[N];
  
  d_constraints = new ArrayList<DistanceConstraint>();
  b_constraints = new  ArrayList<BendingConstraint>();
  
  int hWidth = width/2;
  float dL = (height)/N;
  for(int i=0; i<N; ++i) {    
   // X[i] = new PVector(width/2.0, (i+1)*rest_length);
    X[i] = new PVector( map(i,0,N,hWidth,width), dL);
    F[i] = new PVector(0,0);
    V[i] = new PVector(0,0);
    Ri[i] = new PVector(0,0);
    tmp_X[i] = new PVector(0,0);
    W[i] = (i==0)?0.0f:(1.0f/mass);
    phi0[i]=(float)(Math.PI);
  } 
  
  for(int i=0; i<(N-1); ++i) {    
    AddDistanceConstraint(i, i+1, kStretch);
  }
  
  for(int i=0; i<(N-2); ++i) {    
    AddBendingConstraint(i, i+2, kBend);
  } 
}

public void computeForces() {
   int i=0;
   for(i=0;i<N;i++) {
      F[i].x = 0;
      F[i].y = 0; 

      //add gravity force only for non-fixed points
      if(W[i]!=0)
         F[i].add(gravity); 
   }
}

public void integrateExplicitWithDamping(float deltaTime){
  
  int i=0;

  PVector Xcm = new PVector(0,0);
  PVector Vcm = new PVector(0,0);
  float sumM = 0;
  for(i=0;i<N;i++) {
    V[i] = PVector.mult(V[i], global_damping);                
    V[i] = PVector.add(V[i], PVector.mult(F[i], deltaTime*W[i]));                                            
                
    //calculate the center of mass's position 
    //and velocity for damping calc
    Xcm = PVector.add(Xcm, PVector.mult(X[i],mass));
    Vcm = PVector.add(Vcm, PVector.mult(V[i],mass));
    sumM += mass;
  }
  Xcm=PVector.mult(Xcm, 1.0f/sumM); 
  Vcm=PVector.mult(Vcm, 1.0f/sumM); 

  PMatrix3D I = new PMatrix3D();
  PVector L = new PVector(0,0);
  PVector w = new PVector(0,0);//angular velocity
        
  for(i=0;i<N;i++) { 
     Ri[i] = PVector.sub(X[i], Xcm);   
     L = PVector.add(L, Ri[i].cross(PVector.mult(V[i], mass))); 

     PMatrix3D tmp = new PMatrix3D();
     tmp.set(0.0f,-Ri[i].z,  Ri[i].y, 0.0f,
             Ri[i].z, 0.0f, -Ri[i].x, 0.0f,
             -Ri[i].y, Ri[i].x,0.0f, 0.0f,
             0.0f, 0.0f, 0.0f, 1.0f);
     PMatrix3D tmpT = tmp.get();
     tmpT.transpose();
     tmp.apply(tmpT);
     tmp.scale(mass);
     I.m00 += tmp.m00; I.m01 += tmp.m01; I.m02 += tmp.m02; I.m03 += tmp.m03;
     I.m10 += tmp.m10; I.m11 += tmp.m11; I.m12 += tmp.m12; I.m13 += tmp.m13;
     I.m20 += tmp.m20; I.m21 += tmp.m21; I.m22 += tmp.m22; I.m23 += tmp.m23;
     I.m30 += tmp.m30; I.m31 += tmp.m31; I.m32 += tmp.m32; I.m33 += tmp.m33;
  } 
  I.invert();  
  I.mult(L, w);
        
  //apply center of mass damping
  for(i=0;i<N;i++) {
     PVector delVi = new PVector();
     delVi = PVector.add(Vcm, PVector.sub(w.cross(Ri[i]), V[i]) );               
     V[i]= PVector.add(V[i], PVector.mult(delVi, kDamp));
  }

  //calculate predicted position
  for(i=0;i<N;i++) {
     if(W[i] <= 0.0f) { 
       tmp_X[i] = X[i]; //fixed points
     } else {
       tmp_X[i] = PVector.add(X[i], PVector.mult(V[i],deltaTime));                              
     }
  }  
}

public void stepPhysics() {
  computeForces();
  integrateExplicitWithDamping(dt);
  updateInternalConstraints(dt);
  updateExternalConstraints();
  integrate(dt);
}

public void groundCollision() 
{
  for(int i=0;i<N;i++) {     
     tmp_X[i].x = min(tmp_X[i].x+2, width);
     tmp_X[i].x = max(tmp_X[i].x-2, 0);
     
     tmp_X[i].y = min(tmp_X[i].y+2, height);
     tmp_X[i].y = max(tmp_X[i].y-2, 0);
  }
}


public void updateExternalConstraints() {
   groundCollision();
}


public void updateInternalConstraints(float deltaTime) {
  int i=0;
 
  for (int si=0;si<solver_iterations;++si) {
     for(i=0;i<d_constraints.size();i++) {
         updateDistanceConstraint(i);
     } 
     for(i=0;i<b_constraints.size();i++) {
         updateBendingConstraint(i);
     }         
  }
}

public void updateBendingConstraint(int index) {
  
  BendingConstraint c = b_constraints.get(index); 

  //Using the paper suggested by DevO
  //http://image.diku.dk/kenny/download/kelager.niebe.ea10.pdf
 
  //global_k is a percentage of the global dampening constant 
  float global_k = global_damping*0.01f; 
  PVector center = PVector.mult(PVector.add(tmp_X[c.p1],tmp_X[c.p2]), 0.5f);
  PVector dir_center = PVector.sub(tmp_X[c.p2],center);
  float dist_center = dir_center.mag();

  float diff = 1.0f - ((global_k + c.rest_length) / dist_center);
  PVector dir_force = PVector.mult(dir_center, diff);
  PVector fa = PVector.mult(dir_force, c.k_prime * ((W[c.p1])/c.w));
  PVector fb = PVector.mult(dir_force, c.k_prime * ((W[c.p2])/c.w));
 
  if(W[c.p1] > 0.0f)
     tmp_X[c.p1] = PVector.add(tmp_X[c.p1], fa);
     
  if(W[c.p2] > 0.0f)
     tmp_X[c.p2] = PVector.sub(tmp_X[c.p2], fb);
}

public void updateDistanceConstraint(int t) {
  DistanceConstraint c = d_constraints.get(t);
  PVector dir = PVector.sub(tmp_X[c.p1], tmp_X[c.p2]);
  float len = dir.mag(); 
  if(len <= EPSILON) 
     return;
        
  float w1 = W[c.p1];
  float w2 = W[c.p2];
  float invMass = w1+ w2; 
  if(invMass <= EPSILON) 
     return;
 
  dir.normalize();
  PVector dP = PVector.mult(dir, (1.0f/invMass) * (len-c.rest_length ) * c.k_prime);
  if(w1 > 0.0f)
     tmp_X[c.p1] = PVector.sub(tmp_X[c.p1], PVector.mult(dP, w1));

  if(w2 > 0.0f)
     tmp_X[c.p2] = PVector.add(tmp_X[c.p2], PVector.mult(dP, w2) );   
}

public void integrate(float deltaTime) {
  float inv_dt = 1.0f/deltaTime;
  int i=0; 
  
  for(i=0;i<N;i++) {   
     V[i] = PVector.mult(PVector.sub(tmp_X[i], X[i]), inv_dt);                
     X[i] = tmp_X[i];  
  }
}


public void draw(){
  background(200);
  
  stepPhysics();
  
  //lines
  stroke(127, 120);  
  for(int i=0; i<N-1; ++i) {
    line(X[i].x, X[i].y, X[i+1].x, X[i+1].y);
  }
  
  strokeWeight(6);
  for(int i=0; i<N; ++i) {
    if(intersected==i)
       stroke(0,255,0);
    else
       stroke(0,0,0);
     point(X[i].x, X[i].y);
  } 
}

public void mousePressed() {
  intersected = -1;
  for(int i=0; i<N; ++i) {
     float distX = mouseX-X[i].x;
     float distY = mouseY-X[i].y;
     if( (distX*distX + distY*distY) < 400) {
       intersected = i;
       break;
     }
  }
}

public void mouseDragged() {
   if(intersected != -1) {
      X[intersected].x = mouseX;
      X[intersected].y = mouseY;
      V[intersected].x = 0;
      V[intersected].y = 0;
   }  
}

public void mouseReleased() {
  intersected = -1;
}

public void keyPressed() {
   
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "pbdSpring" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
