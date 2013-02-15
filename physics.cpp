/*

USC/Viterbi/Computer Science
"Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>
#include <utility>
#include <algorithm>
#include <assert.h>

static int collided=0;

// Source: http://rosettacode.org/wiki/Map_range#C.2B.2B
template<typename tVal>
tVal map_value(std::pair<tVal,tVal> a, std::pair<tVal, tVal> b, tVal inVal) {
	tVal inValNorm = inVal - a.first;
	tVal aUpperNorm = a.second - a.first;
	tVal normPosition = inValNorm / aUpperNorm;

	tVal bUpperNorm = b.second - b.first;
	tVal bValNorm = normPosition * bUpperNorm;
	tVal outVal = b.first + bValNorm;

	return outVal;
}

// find value of the force field at an arbitrary position pt inside the bounding box using trilinear interpolation
void computeForceField(struct world* jello, int i, int j, int k, Vec3d* F) {
	point& pt = jello->p[i][j][k];
	{	// map the point to the index of the force field
		std::pair<double,double> box_range(-2.0,2.0), resolution_range(0,jello->resolution-1);
		double i = map_value(box_range, resolution_range, pt.x < -2.0 ? std::max(pt.x, -2.0) : std::min(pt.x, 2.0));
		double j = map_value(box_range, resolution_range, pt.y < -2.0 ? std::max(pt.y, -2.0) : std::min(pt.y, 2.0));
		double k = map_value(box_range, resolution_range, pt.z < -2.0 ? std::max(pt.z, -2.0) : std::min(pt.z, 2.0));
		// lower and upperbound indices of the force fields
		int i0 = ceil(i);
		int i1 = floor(i);
		int j0 = ceil(j);
		int j1 = floor(j);
		int k0 = ceil(k);
		int k1 = floor(k);
		// find weights for each dimension
		double alpha = i - i1, beta = j - j1, gamma = k - k1;
		// assert(0 <= alpha && alpha <= 1.0 && 0 <= beta && beta <= 1.0 && 0 <= gamma && gamma <= 1.0);
		// std::cout << "i:" << i0 << " " << i1 << " " << std::endl;
		// std::cout << "j:" << j0 << " " << j1 << " " << std::endl;
		// std::cout << "k:" << k0 << " " << k1 << " " << std::endl;
		// std::cout << alpha << " " << beta << " " << gamma << std::endl;
		// find force field of the 8 points surrounding the point
		Vec3d* f000 = &jello->forceField[i0 * jello->resolution * jello->resolution + j0 * jello->resolution + k0];
		Vec3d* f001 = &jello->forceField[i0 * jello->resolution * jello->resolution + j0 * jello->resolution + k1];
		Vec3d* f010 = &jello->forceField[i0 * jello->resolution * jello->resolution + j1 * jello->resolution + k0];
		Vec3d* f011 = &jello->forceField[i0 * jello->resolution * jello->resolution + j1 * jello->resolution + k1];
		Vec3d* f100 = &jello->forceField[i1 * jello->resolution * jello->resolution + j0 * jello->resolution + k0];
		Vec3d* f101 = &jello->forceField[i1 * jello->resolution * jello->resolution + j0 * jello->resolution + k1];
		Vec3d* f110 = &jello->forceField[i1 * jello->resolution * jello->resolution + j1 * jello->resolution + k0];
		Vec3d* f111 = &jello->forceField[i1 * jello->resolution * jello->resolution + j1 * jello->resolution + k1];
		// trilinear interpolation
		double fx = 
			(1-alpha)*(1-beta)*(1-gamma)*f000->x +
			(1-alpha)*(beta)*(1-gamma)*f010->x +
			(alpha)*(1-beta)*(1-gamma)*f100->x +
			(alpha)*(beta)*(1-gamma)*f110->x +
			(1-alpha)*(1-beta)*(gamma)*f001->x +
			(1-alpha)*(beta)*(gamma)*f011->x +
			(alpha)*(1-beta)*(gamma)*f101->x +
			(alpha)*(beta)*(gamma)*f111->x;
		double fy = 
			(1-alpha)*(1-beta)*(1-gamma)*f000->y +
			(1-alpha)*(beta)*(1-gamma)*f010->y +
			(alpha)*(1-beta)*(1-gamma)*f100->y +
			(alpha)*(beta)*(1-gamma)*f110->y +
			(1-alpha)*(1-beta)*(gamma)*f001->y +
			(1-alpha)*(beta)*(gamma)*f011->y +
			(alpha)*(1-beta)*(gamma)*f101->y +
			(alpha)*(beta)*(gamma)*f111->y;
		double fz = 
			(1-alpha)*(1-beta)*(1-gamma)*f000->z +
			(1-alpha)*(beta)*(1-gamma)*f010->z +
			(alpha)*(1-beta)*(1-gamma)*f100->z +
			(alpha)*(beta)*(1-gamma)*f110->z +
			(1-alpha)*(1-beta)*(gamma)*f001->z +
			(1-alpha)*(beta)*(gamma)*f011->z +
			(alpha)*(1-beta)*(gamma)*f101->z +
			(alpha)*(beta)*(gamma)*f111->z;
		F->x = fx;
		F->y = fy;
		F->z = fz;
	}
}

void computeCollisionForces(world* jello, int i, int j, int k, point* F) {
	if(jello->p[i][j][k].x > 2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(2.0, jello->p[i][j][k].y, jello->p[i][j][k].z);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
		mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
	if(jello->p[i][j][k].x < -2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(-2.0, jello->p[i][j][k].y, jello->p[i][j][k].z);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
	 	mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
	if(jello->p[i][j][k].z > 2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(jello->p[i][j][k].x, jello->p[i][j][k].y, 2.0);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
		mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
	if(jello->p[i][j][k].z < -2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(jello->p[i][j][k].x, jello->p[i][j][k].y, -2.0);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
	 	mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
	if(jello->p[i][j][k].y > 2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(jello->p[i][j][k].x, 2.0, jello->p[i][j][k].z);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
		mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
	if(jello->p[i][j][k].y < -2.0) {
		/** Compute f_hook **/
		point tmp, L, lNorm;
		point B(jello->p[i][j][k].x, -2.0, jello->p[i][j][k].z);
		mymath::minus(jello->p[i][j][k], B, &L);
		mymath::normalize(L, &lNorm);
		mymath::times(L, -jello->kCollision, &tmp);
	 	mymath::plus(*F, tmp, F);

		/** Compute f_damping **/
		double scalar = -jello->dCollision * mymath::dot(jello->v[i][j][k], lNorm);
		mymath::times(lNorm, scalar, &tmp);
	 	mymath::plus(*F, tmp, F);
	}
}

void computeForces(world* jello, int i, int j, int k, point* F) {
	const double structural_spring_length = 0.142857;
	const double shear_spring_length_type1 = sqrt(2) * 0.142857; // square diagonal shear spring length
	const double shear_spring_length_type2 = sqrt(3) * 0.142857; // cube diagonal shear spring length
	const double bend_spring_length = 2 * 0.142857;
	
	std::vector< Spring > springs;
	
	if(particleExists(i+1, j, k)) {
		springs.push_back(Spring(i+1, j, k, structural_spring_length));
	}
	if(particleExists(i-1, j, k)) {
		springs.push_back(Spring(i-1, j, k, structural_spring_length));
	}
	if(particleExists(i, j+1, k)) {
		springs.push_back(Spring(i, j+1, k, structural_spring_length));
	}
	if(particleExists(i, j-1, k)) {
		springs.push_back(Spring(i, j-1, k, structural_spring_length));
	}
	if(particleExists(i, j, k+1)) {
		springs.push_back(Spring(i, j, k+1, structural_spring_length));
	}
	if(particleExists(i, j, k-1)) {
		springs.push_back(Spring(i, j, k-1, structural_spring_length));
	}

	//Shear: 20 max
	if(particleExists(i, j+1, k+1)) {
		springs.push_back(Spring(i, j+1, k+1, shear_spring_length_type1));
	}
	if(particleExists(i, j+1, k-1)) {
		springs.push_back(Spring(i, j+1, k-1, shear_spring_length_type1));
	}
	if(particleExists(i, j-1, k+1)) {
		springs.push_back(Spring(i, j-1, k+1, shear_spring_length_type1));
	}
	if(particleExists(i, j-1, k-1)) {
		springs.push_back(Spring(i, j-1, k-1, shear_spring_length_type1));
	}
	if(particleExists(i+1, j, k+1)) {
		springs.push_back(Spring(i+1, j, k+1, shear_spring_length_type1));
	}
	if(particleExists(i+1, j, k-1)) {
		springs.push_back(Spring(i+1, j, k-1, shear_spring_length_type1));
	}
	if(particleExists(i-1, j, k+1)) {
		springs.push_back(Spring(i-1, j, k+1, shear_spring_length_type1));
	}
	if(particleExists(i-1, j, k-1)) {
		springs.push_back(Spring(i-1, j, k-1, shear_spring_length_type1));
	}
	if(particleExists(i+1, j+1, k)) {
		springs.push_back(Spring(i+1, j+1, k, shear_spring_length_type1));
	}
	if(particleExists(i+1, j-1, k)) {
		springs.push_back(Spring(i+1, j-1, k, shear_spring_length_type1));
	}
	if(particleExists(i-1, j+1, k)) {
		springs.push_back(Spring(i-1, j+1, k, shear_spring_length_type1));
	}
	if(particleExists(i-1, j-1, k)) {
		springs.push_back(Spring(i-1, j-1, k, shear_spring_length_type1));
	}

	if(particleExists(i+1, j+1, k+1)) {
		springs.push_back(Spring(i+1, j+1, k+1, shear_spring_length_type2));
	}
	if(particleExists(i+1, j+1, k-1)) {
		springs.push_back(Spring(i+1, j+1, k-1, shear_spring_length_type2));
	}
	if(particleExists(i+1, j-1, k+1)) {
		springs.push_back(Spring(i+1, j-1, k+1, shear_spring_length_type2));
	}
	if(particleExists(i+1, j-1, k-1)) {
		springs.push_back(Spring(i+1, j-1, k-1, shear_spring_length_type2));
	}
	if(particleExists(i-1, j+1, k+1)) {
		springs.push_back(Spring(i-1, j+1, k+1, shear_spring_length_type2));
	}
	if(particleExists(i-1, j+1, k-1)) {
		springs.push_back(Spring(i-1, j+1, k-1, shear_spring_length_type2));
	}
	if(particleExists(i-1, j-1, k+1)) {
		springs.push_back(Spring(i-1, j-1, k+1, shear_spring_length_type2));
	}
	if(particleExists(i-1, j-1, k-1)) {
		springs.push_back(Spring(i-1, j-1, k-1, shear_spring_length_type2));
	}

	//Bend: 6 max
	if(particleExists(i+2, j, k)) {
		springs.push_back(Spring(i+2, j, k, bend_spring_length));
	}
	if(particleExists(i-2, j, k)) {
		springs.push_back(Spring(i-2, j, k, bend_spring_length));
	}
	if(particleExists(i, j+2, k)) {
		springs.push_back(Spring(i, j+2, k, bend_spring_length));
	}
	if(particleExists(i, j-2, k)) {
		springs.push_back(Spring(i, j-2, k, bend_spring_length));
	}
	if(particleExists(i, j, k+2)) {
		springs.push_back(Spring(i, j, k+2, bend_spring_length));
	}
	if(particleExists(i, j, k-2)) {
		springs.push_back(Spring(i, j, k-2, bend_spring_length));
	}

	// compute spring forces
    for(std::vector< Spring >::iterator iter = springs.begin(); iter != springs.end(); ++iter) {
		int ni = iter->i, nj = iter->j, nk = iter->k;
		/** Compute f_hook **/
		Vec3d tmp, L, lNorm;
		Vec3d B(jello->p[ni][nj][nk]);
		mymath::minus(jello->p[i][j][k], B, &L); // L = A - B
		double length = mymath::norm(L); // L_bar
		double Rlength;
		Rlength = iter->rest_length; //1/8
		double diff_length = length - Rlength;
		mymath::normalize(L, &lNorm); // L / L_bar
		mymath::times(lNorm, -jello->kElastic * diff_length, &tmp);
		mymath::plus(*F, tmp, F);
		
		/** Compute f_damping **/
		Vec3d va_minus_vb;
		mymath::minus(jello->v[i][j][k], jello->v[ni][nj][nk], &va_minus_vb);
		double scalar = -jello->dElastic * mymath::dot(va_minus_vb, lNorm);
		mymath::times(lNorm, scalar, &tmp);
		mymath::plus(*F, tmp, F);
	}
	
	// compute external force field and collision forces
	Vec3d forceField;
	computeForceField(jello, i, j, k, &forceField);
	Vec3d f_collision;
	computeCollisionForces(jello, i, j, k, &f_collision);
	mymath::plus(forceField, *F, F);
	mymath::plus(f_collision, *F, F);
}

/* Computes acceleration to every control point of the jello cube, 
which is in state given by 'jello'.
Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8]) {
/* for you to implement ... */
	Vec3d f_hook[8][8][8]; // Hook Forces
	Vec3d f_damping[8][8][8]; // Damping Forces
	Vec3d f_forceField[8][8][8]; // Force Field
	
	for(int i = 0; i <= 7; ++i) {
		for(int j = 0; j <= 7; ++j) {
			for(int k = 0; k <= 7; ++k) {
				Vec3d f_total;
				computeForces(jello, i, j, k, &f_total);
				mymath::divide(f_total, jello->mass, &a[i][j][k]);
			}
		}
	}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello) {
	int i,j,k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i=0; i<=7; i++) {
		for (j=0; j<=7; j++) {
			for (k=0; k<=7; k++) {
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
			}
		}
	}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello) {
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];

  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }
  return;  
}
