//============================================================
// STUDENT NAME: Yang Shiyuan
// NUS User ID.: E0518553
// COMMENTS TO GRADER:
//
// ============================================================

#include <cmath>
#include "Sphere.h"

using namespace std;



bool Sphere::hit( const Ray &r, double tmin, double tmax, SurfaceHitRecord &rec ) const 
{
    //***********************************************
    //*********** WRITE YOUR CODE HERE **************
    //***********************************************
    
    // Geometrical Data
    Vector3d localCenter = center;
    double localRadius = radius;

    Vector3d Ro_new = r.origin() - localCenter;
    Vector3d Rd = r.direction();
    Vector3d unit_Rd = Rd.makeUnitVector();

    double a = dot(unit_Rd, unit_Rd);

    //printf("%d\n", a);

    double b = 2.0 * dot(unit_Rd, Ro_new);
    double c = dot(Ro_new, Ro_new) - pow(localRadius, 2.0);

    double d = pow(b, 2.0) - 4 * a * c;

    if (d < 0) 
    {
        //printf("This ray does not touch the surface\n");
        return false;
    }

    double pos_t = (-b + sqrt(d)) / (2.0 * a);
    double neg_t = (-b - sqrt(d)) / (2.0 * a);

    double final_t = 0;
    if (abs(pos_t) >= abs(neg_t)) 
    {
        if (neg_t > 0) 
        {
            final_t = neg_t;
        }
        else 
        {
            final_t = pos_t;
        }
    }

    Vector3d intersect_point = Ro_new + unit_Rd * final_t;
    Vector3d intersect_normal = (intersect_point - localCenter).makeUnitVector();

    if (final_t >= tmin && final_t <= tmax)
    {
        // We have a hit -- populat hit record
        rec.t = final_t;
        rec.p = r.pointAtParam(final_t);
        rec.normal = intersect_normal;
        rec.material = material;
        return true;
    }

    return false;  // YOU CAN REMOVE/CHANGE THIS IF NEEDED.
}



bool Sphere::shadowHit( const Ray &r, double tmin, double tmax ) const 
{
    //***********************************************
    //*********** WRITE YOUR CODE HERE **************
    //***********************************************
    // Geometrical Data
    Vector3d localCenter = center;
    double localRadius = radius;

    Vector3d Ro_new = r.origin() - localCenter;
    Vector3d Rd = r.direction();
    Vector3d unit_Rd = Rd.makeUnitVector();

    double a = dot(unit_Rd, unit_Rd);

    if (a != 1)
    {
        //printf("Check the calculation of a, it should be a unit vector\n");
    }

    double b = 2.0 * dot(unit_Rd, Ro_new);
    double c = dot(Ro_new, Ro_new) - pow(localRadius, 2.0);

    double d = pow(b, 2.0) - 4 * a * c;

    if (d < 0)
    {
        //printf("This ray does not touch the surface\n");
        return false;
    }

    double pos_t = (-b + sqrt(d)) / (2.0 * a);
    double neg_t = (-b - sqrt(d)) / (2.0 * a);

    double final_t = 0;
    if (abs(pos_t) >= abs(neg_t))
    {
        if (neg_t > 0)
        {
            final_t = neg_t;
        }
        else
        {
            final_t = pos_t;
        }
    }

    return (final_t >= tmin && final_t <= tmax);
}
