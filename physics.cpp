/*

USC/Viterbi/Computer Science
"Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>
#include <utility>
#include <assert.h>

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
void findForceField(struct world* jello, const point& pt, Vec3d* F) {
	// map the point to the index of the force field
	std::pair<double,double> box_range(-2.0,2.0), resolution_range(0,jello->resolution-1);
	double i = map_value(box_range, resolution_range, pt.x);
	double j = map_value(box_range, resolution_range, pt.y);
	double k = map_value(box_range, resolution_range, pt.z);
	// lower and upperbound indices of the force fields
	int i0 = ceil(i);
	int i1 = floor(i);
	int j0 = ceil(j);
	int j1 = floor(j);
	int k0 = ceil(k);
	int k1 = floor(k);
	// find weights for each dimension
	double alpha = i - i1, beta = j - j1, gamma = k - k1;
	assert(0 <= alpha && alpha <= 1.0 && 0 <= beta && beta <= 1.0 && 0 <= gamma && gamma <= 1.0);
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
				for(int l = 0; l < jello->p[i][j][k].springs.size(); ++l) {
					/** Compute f_hook **/
					Vec3d tmp, L, lNorm;
					Vec3d B(jello->p[i][j][k].springs.at(l).cp->x,
							jello->p[i][j][k].springs.at(l).cp->y,
							jello->p[i][j][k].springs.at(l).cp->z);
					minus(jello->p[i][j][k], B, &L); // L = A - B
					double length = norm(L); // L_bar
					double rest_length = jello->p[i][j][k].springs.at(l).rest_length; // R
					double diff_length = length - rest_length; // L_bar - R
					normalize(L, &lNorm); // L / L_bar
					times(lNorm, -jello->kElastic * diff_length, &tmp);
					plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);
									
					/** Compute f_damping **/
					Vec3d va_minus_vb;
					minus(jello->v[i][j][k], *(jello->p[i][j][k].springs.at(l).cv), &va_minus_vb);
					double scalar = -jello->dElastic * dot(va_minus_vb, lNorm);
					times(lNorm, scalar, &tmp);
					plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				}
				
				/** Collision **/
				if(jello->p[i][j][k].x > 2.0) {
					/** Compute f_hook **/
					Vec3d tmp, L, lNorm;
					Vec3d B(2.0, jello->p[i][j][k].y, jello->p[i][j][k].z);
					minus(jello->p[i][j][k], B, &L);
					normalize(L, &lNorm);
					times(L, -jello->kCollision, &tmp);
					plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);

					/** Compute f_damping **/
					Vec3d va_minus_vb;
					minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
					double scalar = -jello->dCollision * dot(va_minus_vb, lNorm);
					times(lNorm, scalar, &tmp);
				 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				}
				if(jello->p[i][j][k].x < -2.0) {
					/** Compute f_hook **/
					Vec3d tmp, L, lNorm;
					Vec3d B(-2.0, jello->p[i][j][k].y, jello->p[i][j][k].z);
					minus(jello->p[i][j][k], B, &L);
					normalize(L, &lNorm);
					times(L, -jello->kCollision, &tmp);
					plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);

					/** Compute f_damping **/
					Vec3d va_minus_vb;
					minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
					double scalar = -jello->dCollision * dot(va_minus_vb, lNorm);
					times(lNorm, scalar, &tmp);
				 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				}
				
				// if(jello->p[i][j][k].y > 2.0) {
				// 	Vec3d tmp, L, lNorm;
				// 	Vec3d B(jello->p[i][j][k].x, 2.0, jello->p[i][j][k].z);
				// 	minus(jello->p[i][j][k], B, &L);
				// 	// length of L
				// 	double length = norm(L);
				// 	double rest_length = 0.0;
				// 	double diff_length = length - rest_length;
				// 	normalize(L, &lNorm);
				// 	if(fabs(diff_length) > 0.1) {
				// 		times(lNorm, -jello->kCollision * diff_length, &tmp);
				// 		plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);	
				// 	}
				// 					
				// 	/** Compute f_damping **/
				// 	Vec3d va_minus_vb;
				// 	minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
				// 	double scalar = -jello->dElastic * dot(va_minus_vb, lNorm);
				// 	times(lNorm, scalar, &tmp);
				// 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				// }
				// if(jello->p[i][j][k].y < -2.0) {
				// 	Vec3d tmp, L, lNorm;
				// 	Vec3d B(jello->p[i][j][k].x, -2.0, jello->p[i][j][k].z);
				// 	minus(jello->p[i][j][k], B, &L);
				// 	// length of L
				// 	double length = norm(L);
				// 	double rest_length = 0.0;
				// 	double diff_length = length - rest_length;
				// 	normalize(L, &lNorm);
				// 	if(fabs(diff_length) > 0.1) {
				// 		times(lNorm, -jello->kCollision * diff_length, &tmp);
				// 		plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);	
				// 	}
				// 					
				// 	/** Compute f_damping **/
				// 	Vec3d va_minus_vb;
				// 	minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
				// 	double scalar = -jello->dElastic * dot(va_minus_vb, lNorm);
				// 	times(lNorm, scalar, &tmp);
				// 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				// }
				if(jello->p[i][j][k].z > 2.0) {
					/** Compute f_hook **/
					Vec3d tmp, L, lNorm;
					Vec3d B(jello->p[i][j][k].x, jello->p[i][j][k].y, 2.0);
					minus(jello->p[i][j][k], B, &L);
					normalize(L, &lNorm);
					times(L, -jello->kCollision, &tmp);
					plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);
									
					/** Compute f_damping **/
					Vec3d va_minus_vb;
					minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
					double scalar = -jello->dCollision * dot(va_minus_vb, lNorm);
					times(lNorm, scalar, &tmp);
				 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				}	
				// 	/** Compute f_damping **/
				// 	Vec3d va_minus_vb;
				// 	minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
				// 	double scalar = -jello->dCollision * dot(va_minus_vb, lNorm);
				// 	times(lNorm, scalar, &tmp);
				//  plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				// }
				// if(jello->p[i][j][k].z < -2.0) {
				// 	Vec3d tmp, L, lNorm;
				// 	Vec3d B(jello->p[i][j][k].x, jello->p[i][j][k].y, -2.0);
				// 	minus(jello->p[i][j][k], B, &L);
				// 	// length of L
				// 	double length = norm(L);
				// 	double rest_length = 0.0;
				// 	double diff_length = length - rest_length;
				// 	normalize(L, &lNorm);
				// 	if(fabs(diff_length) > 0.1) {
				// 		times(lNorm, -jello->kCollision * diff_length, &tmp);
				// 		plus(tmp, f_hook[i][j][k], &f_hook[i][j][k]);	
				// 	}
				// 					
				// 	/** Compute f_damping **/
				// 	Vec3d va_minus_vb;
				// 	minus(jello->v[i][j][k], point(0.0, 0.0, 0.0), &va_minus_vb);
				// 	double scalar = -jello->dElastic * dot(va_minus_vb, lNorm);
				// 	times(lNorm, scalar, &tmp);
				// 	plus(tmp, f_damping[i][j][k], &f_damping[i][j][k]);
				// }
				
				/** Computing f_forceField **/
				//findForceField(jello, jello->p[i][j][k], &f_forceField[i][j][k]);
				
				Vec3d f_final;
				plus(f_hook[i][j][k], f_damping[i][j][k], &f_final);
				plus(f_forceField[i][j][k], f_final, &f_final);
				divide(f_final, jello->mass * 512, &a[i][j][k]);
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
