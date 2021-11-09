//============================================================
// STUDENT NAME: Yang Shiyuan
// NUS User ID.: E0518553
// COMMENTS TO GRADER:
//
// ============================================================

#include <cmath>
#include <cfloat>
#include "Vector3d.h"
#include "Color.h"
#include "Ray.h"
#include "Material.h"
#include "Surface.h"
#include "Light.h"
#include "Scene.h"
#include "Raytrace.h"


// This is for avoiding the "epsilon problem" or the shadow acne problem.
static constexpr double DEFAULT_TMIN = 10e-6;

// Use this for tmax for non-shadow ray intersection test.
static constexpr double DEFAULT_TMAX = DBL_MAX;


//////////////////////////////////////////////////////////////////////////////
// Compute the outgoing mirror reflection vector.
// Input incoming vector L is pointing AWAY from surface point.
// Assume normal vector N is unit vector.
// The output reflection vector is pointing AWAY from surface point, and
// has same length as incoming vector L.
//////////////////////////////////////////////////////////////////////////////

static Vector3d mirrorReflect( const Vector3d &L, const Vector3d &N )
{
    return ( 2.0 * dot( N, L ) ) * N - L;
}



//////////////////////////////////////////////////////////////////////////////
// Compute I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vector L, N and V are unit vectors.
//////////////////////////////////////////////////////////////////////////////

static Color computePhongLighting( const Vector3d &L, const Vector3d &N, const Vector3d &V,
                                   const Material &mat, const PointLightSource &ptLight )
{
    Vector3d R = mirrorReflect(L, N);

    auto N_dot_L = (float)dot(N, L);
    if (N_dot_L < 0.0f) N_dot_L = 0.0f;

    auto R_dot_V = (float)dot(R, V);
    if (R_dot_V < 0.0f) R_dot_V = 0.0f;
    float R_dot_V_pow_n = powf(R_dot_V, (float)mat.n);

    return ptLight.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
}



//////////////////////////////////////////////////////////////////////////////
// Traces a ray into the scene.
// reflectLevels: specifies number of levels of reflections (0 for no reflection).
// hasShadow: specifies whether to generate shadows.
//////////////////////////////////////////////////////////////////////////////

Color Raytrace::TraceRay( const Ray &ray, const Scene &scene, 
                          int reflectLevels, bool hasShadow )
{
    Ray uRay( ray );
    uRay.makeUnitDirection();  // Normalize ray direction.


// Find whether and where the ray hits some surface.
// Take the nearest hit point.

    bool hasHitSomething = false;
    double nearest_t = DEFAULT_TMAX;
    SurfaceHitRecord nearestHitRec;

    for (const auto& surface : scene.surfaces)
    {
        SurfaceHitRecord tempHitRec;
        bool hasHit = surface->hit( uRay, DEFAULT_TMIN, DEFAULT_TMAX, tempHitRec );

        if ( hasHit && tempHitRec.t < nearest_t )
        {
            hasHitSomething = true;
            nearest_t = tempHitRec.t;
            nearestHitRec = tempHitRec;
        }
    }

    if ( !hasHitSomething ) return scene.backgroundColor;

    nearestHitRec.normal.makeUnitVector();
    Vector3d N = nearestHitRec.normal;  // Unit vector.
    Vector3d V = -uRay.direction();     // Unit vector.

    Color result( 0.0f, 0.0f, 0.0f );   // The result will be accumulated here.


    //**********************************
    //result = nearestHitRec.material.k_d;  // REMOVE THIS LINE AFTER YOU HAVE FINISHED CODE BELOW.
    //**********************************


// This is to calculate the extremely long term of the Ilocal
// Add to result the phong lighting contributed by each point light source.
// Compute for shadow if hasShadow is true.

    //***********************************************
    //*********** WRITE YOUR CODE HERE **************
    //***********************************************

    Color SumPhongTerm(0.0f, 0.0f, 0.0f);

    for (const auto& ptLight : scene.ptLights) 
    {
        Vector3d unit_L = (ptLight.position - nearestHitRec.p).unitVector();

        Color phongLightResult = computePhongLighting(unit_L, N, V, nearestHitRec.material, ptLight);
        
        if (hasShadow)
        {   
            Ray shadowFeeler(nearestHitRec.p, unit_L);
            bool isBlocked = false;
            Color kShadow(1.0f, 1.0f, 1.0f);
            int surfaceNo = 0;
            for (const auto& surface : scene.surfaces)
            {
                if (surfaceNo < 3)
                {
                    surfaceNo++;
                    continue;
                }

                isBlocked = surface->shadowHit(shadowFeeler, DEFAULT_TMIN, DEFAULT_TMAX);
                if (isBlocked)
                {
                    kShadow.setR(0.0f);
                    kShadow.setG(0.0f);
                    kShadow.setB(0.0f);
                    break;
                }
            }

            phongLightResult *= kShadow;
        }
        
        SumPhongTerm += phongLightResult;
    }
    result += SumPhongTerm;

// This is refering to the calcualtion of Ia * ka -> Local Stuff
// Add to result the global ambient lighting.

    //***********************************************
    //*********** WRITE YOUR CODE HERE **************
    //***********************************************
    Color Ia = scene.amLight.I_a;
    Color Ka = nearestHitRec.material.k_a;
    Color IaxKa = Ia*Ka;

    result += IaxKa;


// This is refering to the calculation of Krg * Ireflected
// Add to result the reflection of the scene.

    //***********************************************
    //*********** WRITE YOUR CODE HERE **************
    //***********************************************

    if (reflectLevels == 0) 
    {
        return result;
    } 

    Color Krg = nearestHitRec.material.k_rg;
    Vector3d R = mirrorReflect(V, N);
    Ray reflectRay(nearestHitRec.p , R);
    return result += Krg * TraceRay(reflectRay, scene, reflectLevels-1, hasShadow);
}
