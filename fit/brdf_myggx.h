#ifndef _BRDF_MYGGX_
#define _BRDF_MYGGX_

#include "brdf.h"

class BrdfMyGGX : public Brdf
{
public:
    virtual float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
    {
        if (V.z <= 0)
        {
            pdf = 0;
            return 0;
        }

        // // masking
        // const float LambdaV = lambda(alpha, V.z);

        // shadowing
        float G2;
        if (L.z <= 0.0f)
            G2 = 0;
        else
        {
            // const float LambdaL = lambda(alpha, L.z);
            // G2 = 1.0f/(1.0f + LambdaV + LambdaL);

            G2 = lambda(alpha, L.z) * lambda(alpha, V.z);
        }

        // D
        // This is GGX, but it terms out the alpha=roughness^2 here is the roughness of our implementation
        // So we modify the alpha^2 here to single alpha
        /*
D = alpha^2 / (alpha^2. + (H.x^2 /H.z)^2 + (H.y/H.z)^2))
  = alpha^2 H.z^2 / (alph^2H.z^2 + H.x^2 + H.y^2)
  = alpha^2 H.z^2 / ((alpha^2 - 1) H.z^2 + 1)
D^2 = alpha^4 H.z^4 / ((alpha^2 - 1) H.z^2 + 1)^2
D^2 / (pi * alpha^2 H.z^4) = alpha^2 / pi /  (alpha^2 - 1) H.z^2 + 1) ^ 2
        */

        const vec3 H = normalize(V + L);

        // const float slopex = H.x/H.z;
        // const float slopey = H.y/H.z;
        // // float D = 1.0f / (1.0f + (slopex*slopex + slopey*slopey)/alpha/alpha);
        // float D = 1.0f / (1.0f + (slopex*slopex + slopey*slopey)/alpha);
        // D = D*D;
        // // D = D/(3.14159f * alpha*alpha * H.z*H.z*H.z*H.z);
        // D = D/(3.14159f * alpha * H.z*H.z*H.z*H.z);

        float NdotH = std::max(0.f, H.z);
        float d = ((NdotH * alpha - NdotH) * NdotH + 1);
        float D = alpha / (d * d * 3.14159265359);

        float VdotH = std::max(0.f, dot(V, H));
        float inv_VdotH = VdotH <= 1e-6 ? 0.0 : (1.0 / VdotH);
        float inv_Vz = V.z <= 1e-6 ? 0.0 : (1.0 / V.z);

        pdf = fabsf(D * H.z / 4.0f * inv_VdotH);
        float res = D * G2 / 4.0f * inv_Vz;

        return res;
    }

    virtual vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
    {
        // const float phi = 2.0f*3.14159f * U1;
        // const float r = alpha*sqrtf(U2/(1.0f - U2));
        // const vec3 N = normalize(vec3(r*cosf(phi), r*sinf(phi), 1.0f));
        // const vec3 L = -V + 2.0f * N * dot(N, V);
        // return L;
        const float phi = 2.0f*3.14159f * U1;

        const float cosThetaH = sqrtf(std::max(0.0f, (1.0f-U2)/((alpha-1.0f)*U2+1.0f) ));
        const float sinThetaH = sqrtf(std::max(0.0f, 1.0f - cosThetaH * cosThetaH));

        // Get our GGX NDF sample (i.e., the half vector)
        vec3 H = vec3(sinThetaH * cos(phi), (sinThetaH * sin(phi)), cosThetaH);
        vec3 L = 2 * dot(V, H) * H - V;
        return L;
    }

private:
    float lambda(const float alpha, const float cosTheta) const
    {
        // const float a = 1.0f / alpha / tanf(acosf(cosTheta));
        // return (cosTheta < 1.0f) ? 0.5f * (-1.0f + sqrtf(1.0f + 1.0f/a/a)) : 0.0f;    
        
        // Previously this is TR, but we use simpler Schlick
        float k = alpha / 2.0; // Simplified Schlick
        return cosTheta / (cosTheta * (1.0f - k) + k);
    }
};

#endif
