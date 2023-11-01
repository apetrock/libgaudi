#ifndef __GGSHADERS__
#define __GGSHADERS__

#include <iostream>
#include <string>
namespace gg {
inline std::string get_shader(std::string name, int width = 1280,
                              int height = 720) {
  std::string f_buff_vert =
      R"(
#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec2 aTexCoords;
layout(location = 2) in vec3 aNormal;
layout(location = 3) in vec3 aColor;

// Output data ; will be interpolated for each fragment.
out vec2 UV;
out vec3 Position_worldspace;
out vec3 Normal_cameraspace;
out vec3 Color_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;

// Values that stay constant for the whole mesh.
uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
uniform vec3 LightPosition_worldspace;

out vec3 ec_pos;

void main(){

  //use for the normals
  mat4 MVP = P * V * M;
  ec_pos = (M * V * vec4(aPos,1)).xyz;	
  // Output position of the vertex, in clip space : MVP * position
  gl_Position =  MVP * vec4(aPos,1);
	
  // Position of the vertex, in worldspace : M * position
  Position_worldspace = (M * vec4(aPos,1)).xyz;
	
  // Vector that goes from the vertex to the camera, in camera space.
  // In camera space, the camera is at the origin (0,0,0).
  vec3 vertexPosition_cameraspace = ( V * M * vec4(aPos,1)).xyz;
  EyeDirection_cameraspace = vec3(0,0,0) - vertexPosition_cameraspace;

  // Vector that goes from the vertex to the light, in camera space. M is ommited because it's identity.
  vec3 LightPosition_cameraspace = ( V * vec4(LightPosition_worldspace,1)).xyz;
  LightDirection_cameraspace = LightPosition_cameraspace + EyeDirection_cameraspace;
	
  // Normal of the the vertex, in camera space
  Normal_cameraspace = ( V * M * vec4(aNormal,0)).xyz; // Only correct if ModelMatrix does not scale the model ! Use its inverse transpose if not.
	Color_cameraspace = aColor;
  // UV of the vertex. No special space for this one.
  UV = aTexCoords;
}
)";

  std::string f_buff_frag =
      R"(

#version 330 core

// Interpolated values from the vertex shaders
in vec2 UV;
in vec3 Position_worldspace;
in vec3 Normal_cameraspace;
in vec3 Color_cameraspace;
in vec3 EyeDirection_cameraspace;
in vec3 LightDirection_cameraspace;

// Ouput data
out vec4 color;

in vec3 ec_pos;
// Values that stay constant for the whole mesh.
uniform vec3 LightPosition_worldspace;
void main(){

	// Light emission properties
	// You probably want to put them as uniforms
	vec3 LightColor = vec3(1,1,1);
	float LightPower = 30.0f;
	
	// Material properties
	vec3 MaterialDiffuseColor = Color_cameraspace;
	vec3 MaterialAmbientColor = vec3(0.25,0.25,0.25) * MaterialDiffuseColor;
	vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);

	// Distance to the light
	float distance = length( LightPosition_worldspace - Position_worldspace );

	// Normal of the computed fragment, in camera space
	
	vec3 ec_normal = normalize(cross(dFdx(ec_pos), dFdy(ec_pos)));
	vec3 n = ec_normal;
	//vec3 n = normalize( Normal_cameraspace );
	// Direction of the light (from the fragment to the light)
	vec3 l = normalize( LightDirection_cameraspace );
	// Cosine of the angle between the normal and the light direction, 
	// clamped above 0
	//  - light is at the vertical of the triangle -> 1
	//  - light is perpendicular to the triangle -> 0
	//  - light is behind the triangle -> 0
	float cosTheta = clamp( dot( n,l ), 0,1 );
	
	// Eye vector (towards the camera)
	vec3 E = normalize(EyeDirection_cameraspace);
	// Direction in which the triangle reflects the light
	vec3 R = reflect(-l,n);
	// Cosine of the angle between the Eye vector and the Reflect vector,
	// clamped to 0
	//  - Looking into the reflection -> 1
	//  - Looking elsewhere -> < 1
	float cosAlpha = clamp( dot( E,R ), 0,1 );
	
	//color = vec3(cosTheta);
		//Ambient : simulates indirect lighting
	vec3 color3 = MaterialAmbientColor +
	  //Diffuse : "color" of the object
	  MaterialDiffuseColor * LightColor * LightPower * cosTheta/(distance*distance) + 
	
      	// Specular : reflective highlight, like a mirror
	  MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distance*distance);
  //color = vec4(UV, 0.0,1.0);
  color = vec4(color3,1.0);
}
)"; //

  //////////////////////////////////////////////////////
  std::string g_buff_vert = //
      R"(
  #version 330 core
  layout(location = 0) in vec3 aPos;
  layout(location = 1) in vec2 aTexCoords;
  layout(location = 2) in vec3 aNormal;
  layout(location = 3) in vec3 aColor;

  out vec3 FragPos;
  out vec2 TexCoords;
  out vec3 Normal;
  out vec3 Color;
  
  uniform mat4 P;
  uniform mat4 V;
  uniform mat4 M;

  void main()
  {   
      mat4 MVP = P * V * M;

      vec4 worldPos = V * M * vec4(aPos, 1.0);
      FragPos = worldPos.xyz; 
      TexCoords = aTexCoords;
      
      mat3 normalMatrix = transpose(inverse(mat3(M)));
      Normal = normalMatrix * aNormal;
      Color = aColor;
      //gl_Position =  MVP * vec4(aPos,1);
      gl_Position = P * worldPos;
  }
)";

  std::string g_buff_frag =
      R"(
  #version 330 core
  layout (location = 0) out vec3 gPosition;
  layout (location = 1) out vec3 gNormal;
  layout (location = 2) out vec4 gAlbedoSpec;

  in vec2 TexCoords;
  in vec3 FragPos;
  in vec3 Normal;
  in vec3 Color;
  
  uniform sampler2D texture_diffuse1;
  uniform sampler2D texture_specular1;

  void main()
  {    
      // store the fragment position vector in the first gbuffer texture
      gPosition = FragPos;
      // also store the per-fragment normals into the gbuffer
      gNormal = normalize(Normal);
      gNormal = normalize(cross(dFdx(FragPos), dFdy(FragPos)));
	
      // and the diffuse per-fragment color
      gAlbedoSpec.rgb = texture(texture_diffuse1, TexCoords).rgb;
      // store specular intensity in gAlbedoSpec's alpha component
      gAlbedoSpec.a = texture(texture_specular1, TexCoords).r;
      gAlbedoSpec = vec4(Color,1.0);
  }
)";

  //////////////////////////////////////////////////////
  std::string sqr_vert =
      R"(
#version 330 core
layout (location = 0) in vec3 aPos;

out vec2 TexCoords;

void main()
{
    TexCoords = aPos.xy * 0.5 - 0.5;
    gl_Position = vec4(aPos, 1.0);
}
)";

  std::string hdr_sqr_frag =
      R"(
#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D hdrBuffer;
uniform bool hdr;
uniform float exposure;

void main()
{   
    const float gamma = 2.2;
    vec3 hdrColor = texture(hdrBuffer, TexCoords).rgb;
    //FragColor = vec4(hdrColor, 1.0);
    //return;

    if(hdr)
    {
        // reinhard
        // vec3 result = hdrColor / (hdrColor + vec3(1.0));
        // exposure
        vec3 result = vec3(1.0) - exp(-hdrColor * exposure);
        // also gamma correct while we're at it       
        result = pow(result, vec3(1.0 / gamma));
        FragColor = vec4(result, 1.0);
    }
    else
    {
        vec3 result = pow(hdrColor, vec3(1.0 / gamma));
        FragColor = vec4(result, 1.0);
    }
}
)";

  std::string def_sqr_frag =
      R"(
        #version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D gAlbedoSpec;
uniform sampler2D ssao;
uniform sampler2D bleed;

uniform vec3 viewPos;

void main()
{             
    // retrieve data from gbuffer

    vec3 FragPos = texture(gPosition, TexCoords).rgb;
    vec3 Normal = texture(gNormal, TexCoords).rgb;
    vec3 Diffuse = texture(gAlbedoSpec, TexCoords).rgb;
    float AmbientOcclusion = texture(ssao, TexCoords).r;
    vec3 Bleed = texture(bleed, TexCoords).rgb;
    
    //FragColor = vec4(Normal, 1.0);
    //FragColor = vec4(vec3(AmbientOcclusion), 1.0);
    //FragColor = vec4(Bleed, 1.0);
    //return;
    
    float Specular = 0.5;
    // then calculate lighting as usual

    vec3 ambient  =  Diffuse; // hard-coded ambient component
    vec3 lighting  = vec3(0.0); // hard-coded ambient component
    //FragColor = vec4(lighting, 1.0);
    //return;

    vec3 viewDir  = normalize(viewPos - FragPos);
    //for(int i = 0; i < NR_LIGHTS; ++i)
    for(int i = 0; i < 1; ++i)
    
    {
        // diffuse
        //vec3 lpos = lights[i].Position;
        //float llin = lights[i].Linear;
        //float lcol = lights[i].Color;
        
        vec3 lcol = 5.0 * vec3(1.0, 1.0, 1.0);
        vec3 lpos = vec3(1.0, 1.0, 1.0);
        float llin = 0.05;
        float lquad = 0.8;

        vec3 lightDir = normalize(lpos - FragPos);


        vec3 diffuse = max(dot(Normal, lightDir), 0.0) * Diffuse * lcol;
        // specular
        vec3 halfwayDir = normalize(lightDir + viewDir);  
        float spec = pow(max(dot(Normal, halfwayDir), 0.0), 16.0);

        vec3 specular = lcol * spec * Specular;
        
        // attenuation
        float distance = length(lpos - FragPos);
        float attenuation = 1.0 / (1.0 + llin * distance + lquad * distance * distance);
        diffuse *= attenuation;
        specular *= attenuation;
        lighting += diffuse + specular;
    }


    vec2 tc = 1.0 + TexCoords;
    //if(tc.x < 0.5)
      FragColor = vec4((1.0 * lighting + 0.8 * AmbientOcclusion * (1.5 * Bleed + 0.2 * ambient)), 1.0);
    //else
    //  FragColor = vec4(vec3(8.0 * Bleed), 1.0);
}
)";

  std::string ssao_sqr_frag =
      R"(
        #version 330 core

#define KSIZE 256
#define PSI  1.533751168755204288118041 
#define PHI  sqrt(2.0)
#define ROOT3  sqrt(3.0)

#define GPHI 0.5 * (1.0 +  sqrt(5.0))

#define LIGHT_POWER 20.
#define SPECULAR_POWER 20.
#define AMBIENT .3

#define PI     3.14159265
#define TWO_PI 6.28318530

out float FragColor;

in vec2 TexCoords;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
//uniform sampler2D texNoise;

uniform vec3 samples[256];

// parameters (you'd probably want to use them as uniforms to more easily tweak the effect)

float radius = 0.2;
float bias = -0.01;

// tile noise texture over screen based on screen dimensions divided by noise size
const vec2 noiseScale = vec2(float(SCREEN_WIDTH)/4.0, float(SCREEN_HEIGHT)/4.0); 

struct quat
{
    float s;
    vec3 v;
};

vec3 rotate(quat q, vec3 p) // NOTE: order of parameters copies order of applying rotation matrix: M v
{
    return p + 2.0 * cross(q.v, cross(q.v, p) + q.s * p); // suggested by mla, requires q to be unit (i.e. normalized)

    // Derive to https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula#Vector_formulation
    //return p + 2.0 * cross(q.v, cross(q.v, p)) + 2.0 * cross(q.v, q.s * p); // cross-product is distributive
    //return p + 2.0 * cross(q.v, q.s * p) + 2.0 * cross(q.v, cross(q.v, p)); // vector addition is commutative
    //return p + 2.0 * q.s * cross(q.v, p) + 2.0 * cross(q.v, cross(q.v, p)); // scalar can be factored-out
    // translate variable names
    vec3 x = p;
    float a = q.s;
    vec3 omega = q.v;
    return x + 2.0 * a * cross(omega, x) + 2.0 * cross(omega, cross(omega, x)); // Euler Rodrigues' Formula
}

uniform mat4 projection;
vec3 getSphere(int i){
    float N = float(KSIZE);
    float t = float(i) / N;
    float phi = GPHI;
    float psi = ROOT3;
    float sqti = sqrt(t);
    float sqt1i = sqrt(1.0 - t);
    float thet = 2.0 * PI *N * t;
    float NtP = (N * t) / PSI;
    float t0 = NtP - floor(NtP);
    float x0 =  sqti * sin(thet / PHI);
    float x1 =  sqti * cos(thet / PHI);
    float y0 =  sqt1i * sin(thet / PSI);
    float y1 =  sqt1i * cos(thet / PSI);
    quat q = quat(y0, vec3(y1, x0, x1));
    vec3 p0 = 1.0 * pow(t0, 2.5/3.0) * vec3(1.0, 0.0, 0.0);
    return rotate(q, p0);
}

float gold_noise(in vec2 xy, in float seed)
{
  return fract(tan(distance(xy*GPHI, xy)*seed)*xy.x);
}

vec3 rvec(vec2 t){
  return vec3 (gold_noise(t, 0.0 * GPHI),  // r
                    gold_noise(t, 1.0 * GPHI),  // g
                    gold_noise(t, 2.0* GPHI)); // b
}

void main()
{

    // get input for SSAO algorithm
    vec3 fragPos = texture(gPosition, TexCoords).xyz;
    vec3 normal = normalize(texture(gNormal, TexCoords).rgb);
    vec3 randomVec = 0.5 - normalize(rvec(500.0 * TexCoords));
    //vec3 randomVec = normalize(texture(texNoise, TexCoords * noiseScale).xyz);
    
    // create TBN change-of-basis matrix: from tangent-space to view-space
    vec3 tangent = normalize(randomVec - normal * dot(randomVec, normal));
    vec3 bitangent = cross(normal, tangent);

    mat3 TBN = mat3(tangent, bitangent, normal);
    // iterate over the sample kernel and calculate occlusion factor
    float occlusion = 0.0;

    for(int i = 0; i < KSIZE; ++i)
    {
        // get sample position
        //vec3 samplePos = TBN * samples[i]; // from tangent to view-space
        vec3 sphere = getSphere(i);
        sphere.z = abs(sphere.z);
        vec3 samplePos = TBN * sphere; // from tangent to view-space
        
        samplePos = fragPos + samplePos * radius; 
        
        // project sample position (to sample texture) (to get position on screen/texture)
        vec4 offset = vec4(samplePos, 1.0);
        offset = projection * offset; // from view to clip-space
        offset.xyz /= offset.w; // perspective divide
        offset.xyz = offset.xyz * 0.5 + 0.5; // transform to range 0.0 - 1.0

        //occlusion += texture(gPosition, offset.xy).z;

        // get sample depth
        float sampleDepth = texture(gPosition, offset.xy).z; // get depth value of kernel sample
        
        // range check & accumulate
        float rangeCheck = smoothstep(0.0, 1.0, radius / abs(fragPos.z - sampleDepth));
        //occlusion += (sampleDepth >= samplePos.z + bias ? 1.0 : 0.0);           
        //occlusion += rangeCheck;
        occlusion += (sampleDepth >= samplePos.z + bias ? 1.0 : 0.0) * rangeCheck;           
    }
    //occlusion /= kernelSize;
    occlusion = 1.0 - (occlusion / float(KSIZE));
    
    FragColor = pow(occlusion, 1.8);
}
)";

  std::string bleed_sqr_frag =
      R"(
        #version 330 core

#define KSIZE 128
#define PSI  1.533751168755204288118041 
#define PHI  sqrt(2.0)
#define ROOT3  sqrt(3.0)

#define GPHI 0.5 * (1.0 +  sqrt(5.0))

#define LIGHT_POWER 20.
#define SPECULAR_POWER 20.
#define AMBIENT .3

#define PI     3.14159265
#define TWO_PI 6.28318530

out vec4 FragColor;

in vec2 TexCoords;

uniform mat4 projection;
uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D gAlbedoSpec;
//uniform sampler2D texNoise;
uniform vec3 samples[256];
uniform vec3 viewPos;

// parameters (you'd probably want to use them as uniforms to more easily tweak the effect)

float radius = 0.1;
float bias = 0.02;

// tile noise texture over screen based on screen dimensions divided by noise size
const vec2 noiseScale = vec2(float(SCREEN_WIDTH)/4.0, float(SCREEN_HEIGHT)/4.0); 

struct quat
{
    float s;
    vec3 v;
};
/*
vec3 getSphere(int i){
    float N = float(KSIZE);
    float t = float(i) / N;
    float sqti = sqrt(t);
    float phi = GPHI;
    float Nt = N * t;
    float thet = 2.0 * PI * Nt;
        
    float x0 = sqti * cos(thet * phi);
    float y0 = sqti * sin(thet * phi);
    return vec3(x0, y0, 1.0);
}
*/
vec3 rotate(quat q, vec3 p) // NOTE: order of parameters copies order of applying rotation matrix: M v
{
    return p + 2.0 * cross(q.v, cross(q.v, p) + q.s * p); // suggested by mla, requires q to be unit (i.e. normalized)
}

vec3 getSphere(int i){
    float N = float(KSIZE);
    float t = float(i) / N;
    float phi = GPHI;
    float psi = ROOT3;
    float sqti = sqrt(t);
    float sqt1i = sqrt(1.0 - t);
    float thet = 2.0 * PI *N * t;
    float NtP = (N * t) / PSI;
    float t0 = NtP - floor(NtP);
    float x0 =  sqti * sin(thet / PHI);
    float x1 =  sqti * cos(thet / PHI);
    float y0 =  sqt1i * sin(thet / PSI);
    float y1 =  sqt1i * cos(thet / PSI);
    quat q = quat(y0, vec3(y1, x0, x1));
    vec3 p0 = 1.0 * pow(t0, 2.5/3.0) * vec3(1.0, 0.0, 0.0);
    return rotate(q, p0);
}

float gold_noise(in vec2 xy, in float seed)
{
  return fract(tan(distance(xy*GPHI, xy)*seed)*xy.x);
}

vec3 rvec(vec2 t){
  return vec3 (gold_noise(t, 0.0 * GPHI),  // r
                    gold_noise(t, 1.0 * GPHI),  // g
                    gold_noise(t, 2.0* GPHI)); // b
}

vec3 reflect(const vec3 N, vec3 x) {
  vec3 rx = x - 2.0 * dot(x, N) * N;
  return rx;
};

vec3 project(vec3 pos){
    vec4 offset = vec4(pos, 1.0);
    offset = projection * offset; // from view to clip-space
    offset.xyz /= offset.w; // perspective divide
    offset.xyz = offset.xyz * 0.5 + 0.5; // transform to range 0.0 - 1.0
    return offset.xyz;
}

void find_pos(vec3 p0, vec3 p1, out vec3 gpos, out vec3 gnorm){
    vec3 dir = normalize(p1 - p0);
    vec3 offset = project(p1);

    gpos = texture(gPosition, offset.xy).rgb;
    gnorm = texture(gNormal, offset.xy).rgb;

    vec3 cdir = gpos - p0;
    float d = dot(cdir, dir);
    vec3 p1p = p0 + d * normalize(cdir);
    offset = project(p1p);

    gpos = texture(gPosition, offset.xy).rgb;
    gnorm = texture(gNormal, offset.xy).rgb;

}

void main()
{

    // get input for SSAO algorithm
    vec3 fragPos = texture(gPosition, TexCoords).xyz;
    vec3 normal = normalize(texture(gNormal, TexCoords).rgb);
    //vec3 normal = vec3(0.0, 0.0, 1.0);
    vec3 randomVec = normalize(rvec(499.0 * TexCoords));
    //vec3 randomVec = normalize(texture(texNoise, TexCoords * noiseScale).xyz);
    
    // create TBN change-of-basis matrix: from tangent-space to view-space
    vec3 tangent = normalize(randomVec - normal * dot(randomVec, normal));
    vec3 bitangent = cross(normal, tangent);

    mat3 TBN = mat3(tangent, bitangent, normal);
    // iterate over the sample kernel and calculate occlusion factor
    vec3 occlusion = vec3(0.0);
    vec3 viewDir  = normalize(viewPos - fragPos);
    float accum = 0.0;
    for(int i = 0; i < KSIZE; ++i)
    {
        // get sample position
        //vec3 samplePos = TBN * samples[i]; // from tangent to view-space
        vec3 sphere = getSphere(i);
        sphere.z = abs(sphere.z);
        vec3 samplePos = TBN * sphere; // from tangent to view-space

        samplePos = fragPos + samplePos * radius; 

        vec3 p0, p1, p2;
        float cdist;
        //guess
        p0 = viewPos;
        p1 = samplePos;
        
        vec3 gpos, gnorm0, gnorm1, offset;
        
        offset = project(p1);
        gpos = texture(gPosition, offset.xy).rgb;
        gnorm0 = texture(gNormal, offset.xy).rgb;
        float sampleDepth = gpos.z; // get depth value of kernel sample
        float rc = smoothstep(0.0, 1.0, radius / abs(fragPos.z - sampleDepth));

        if(gpos.z <= p1.z + bias){

          vec3 dv = -normalize(reflect(gnorm0, viewDir));
          float c0 = dot(viewDir, gnorm0);
          occlusion += c0 * texture(gAlbedoSpec, offset.xy).rgb;
          continue;
        } 
        continue;
        //find_pos(p0, p1, gpos, gnorm0);
        p1 = gpos;

        offset = project(p1);

        vec3 d0 = normalize(p1 - p0);

        //calc new direction
        vec3 d1 = -normalize(reflect(gnorm0, d0));
        p2 = p1 + (0.25 * radius  + 0.75 * gold_noise(250.0 * TexCoords, GPHI) * radius) * d1;

        offset = project(p2);

        gpos = texture(gPosition, offset.xy).rgb;
        gnorm1 = texture(gNormal, offset.xy).rgb;
        if(gpos.z < p2.z + bias) continue;

        find_pos(p1, p2, gpos, gnorm1);

        float c0 = -dot(d1, gnorm0);
        float c1 = -dot(d1, gnorm1);
        c0 = clamp(c0, 0.0, 1.0);
        c1 = clamp(c1, 0.0, 1.0);

        p2 = gpos;
        float dist = length(p2 - p1); //use dot prod instea

        float denom = 1.0 / (dist * dist + 1e-6);
        denom = clamp(denom, 0.0, 1.0);
        
        occlusion += c0 * c1 * denom * texture(gAlbedoSpec, offset.xy).rgb;
    }

    float A = PI * radius * radius / float(KSIZE);
    //occlusion *= A;
    occlusion /= float(KSIZE);
    
    FragColor = vec4(1.0 * occlusion, 1.0);
}
)";

  std::string blur_sqr_frag =
      R"(
        #version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D targetInput;
uniform sampler2D gNormal;

void main() 
{
    vec2 texelSize = 1.0 / vec2(textureSize(targetInput, 0));
    vec3 dx = vec3(0.0);
    float accum = 0.0;
    vec3 N0 = texture(gNormal, TexCoords).xyz;
    for (int x = -2; x < 2; ++x) 
    {
        for (int y = -2; y < 2; ++y) 
        {
            vec2 offset = vec2(float(x), float(y)) * texelSize;
            vec3 Ni = texture(gNormal, TexCoords + offset).xyz;
            float dNN = dot(N0, Ni) + 5e-2;

            accum += 1.0;
            dx += dNN * ( texture(targetInput, TexCoords + offset)).xyz;
        }
    }
    
    FragColor = vec4(dx, 1.0)/accum;
}  
)";

  //////////////////////////////////////////////////////
  std::string dbg_point_vert =
      R"(
        #version 330 core
        layout(location = 0) in vec3 vertexPosition_modelspace;
        layout(location = 1) in vec2 vertexUV;
        layout(location = 3) in vec3 vertexColor_modelspace;
        out vec2 UV;
        out vec3 Color_cameraspace;
        uniform mat4 MVP;

        void main() {

            Color_cameraspace = vertexColor_modelspace;
            UV = vertexUV;
            gl_Position = MVP * vec4(vertexPosition_modelspace, 1.0);
        }
)";

  std::string dbg_point_frag =
      R"(
        /* Fragment shader */
        #version 330
        in vec2 UV;
        in vec3 Position_worldspace;
        in vec3 Color_cameraspace;

        out vec4 color;
        void main() {
            color = vec4(Color_cameraspace, 1.0);
        }
)";

  std::string dbg_point_geo =
      R"(

)";

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  std::string dbg_line_vert =
      R"(
        #version 330 core
        layout(location = 0) in vec3 aPos;
        layout(location = 1) in vec3 aColor;

        out vec3 color_vert;
        out vec3 wpos_vert;
        uniform mat4 P;
        uniform mat4 M;
        uniform mat4 V;

        void main() {
            mat4 MVP = P * V * M;
            color_vert = aColor;
            vec4 wpos = V * M * vec4(aPos, 1.0);
            wpos_vert = wpos.xyz;
            gl_Position = P * wpos;
        }
)";

  std::string dbg_line_geo =
      R"(
/* geometry shader */
#version 330 core
layout (lines) in;
layout (triangle_strip, max_vertices = 128) out;

uniform mat4 P;
uniform mat4 M;
uniform mat4 V;

in vec3 color_vert[];
in vec3 wpos_vert[];

out vec3 wpos_frag;
out vec3 color_frag;
void main() {
  mat4 MVP = P * V * M;

  float r = 0.003;
    vec3 p0 = wpos_vert[0];
    vec3 p1 = wpos_vert[1];
    
    vec3 T = normalize(p1 - p0);
    vec3 c = 0.5 * (p1 + p0);
    //vec3 N = -normalize(c);

    vec3 B = normalize(cross(T, normalize(c)));
    vec3 N = normalize(cross(B, T));
    
    const int Nv = 6;
    float Nvf = float(Nv);
    vec3[Nv] circ = vec3[Nv]( 0.0 );
    for (int i = 0; i < Nv; i++) {
        float thet = float(i) / Nvf * 2.0 * 3.14159265359;
        float s0 = sin(thet);
        float c0 = cos(thet);
        circ[i] = s0 * N + c0 * B;    
    }

    //endcaps
    
    const int Nk = 3;
    const float Nkf = float(Nk);
    for(int k = 0; k < Nk; k++){
      float z0 = float(k) * r / Nkf;
      float z1 = float(k+1) * r / Nkf;
      for (int i0 = Nv-1; i0 > -1; i0--) {
            int i1 = i0 % Nv;
            float r0 = sqrt(r*r -  z0*z0);
            float r1 = sqrt(r*r -  z1*z1);
            vec3 t0 = p0 - (z0 * T + r0 * circ[i1]);
            vec3 t1 = p0 - (z1 * T + r1 * circ[i1]);
            
            wpos_frag = t0;
            gl_Position = P * vec4(t0, 1.0);
            color_frag = color_vert[0];
            EmitVertex();

            wpos_frag = t1;
            gl_Position = P * vec4(t1, 1.0);
            color_frag = color_vert[1];
            EmitVertex();
          }
    }
    
    for(int k = 0; k < Nk; k++){
      float z0 = float(k) * r / Nkf;
      float z1 = float(k+1) * r / Nkf;
      for (int i0 = 0; i0 < Nv + 1; i0++) {
            int i1 = i0 % Nv;
            float r0 = sqrt(r*r -  z0*z0);
            float r1 = sqrt(r*r -  z1*z1);
            vec3 t0 = p1 + z0 * T - r0 * circ[i1];
            vec3 t1 = p1 + z1 * T - r1 * circ[i1];
            
            wpos_frag = t0;
            gl_Position = P * vec4(t0, 1.0);
            color_frag = color_vert[0];
            EmitVertex();

            wpos_frag = t1;
            gl_Position = P * vec4(t1, 1.0);
            color_frag = color_vert[1];
            EmitVertex();
          }
    }

    for (int i0 = 0; i0 < Nv + 1; i0++) {
        int i1 = i0 % Nv;
        vec3 t0 = p0 + r * circ[i1];
        vec3 t1 = p1 + r * circ[i1];
        
        wpos_frag = t0;
        gl_Position = P * vec4(t0, 1.0);
        color_frag = color_vert[0];
        EmitVertex();

        wpos_frag = t1;
        gl_Position = P * vec4(t1, 1.0);
        color_frag = color_vert[1];
        EmitVertex();
    }

    EndPrimitive();
}
)";

  std::string dbg_line_frag =
      R"(
        /* Fragment shader */
        #version 330              
        layout (location = 0) out vec3 gPosition;
        layout (location = 1) out vec3 gNormal;
        layout (location = 2) out vec4 gAlbedoSpec;

        in vec2 UV;
        in vec3 Position_worldspace;
        in vec3 color_frag;
        in vec3 wpos_frag;

        void main() {
            
          	vec3 ec_normal = normalize(cross(dFdx(wpos_frag), dFdy(wpos_frag)));
            gAlbedoSpec = vec4(color_frag, 1.0);
            gNormal = ec_normal;
            gPosition = wpos_frag;
        }
)";

  if (name == "f_buff_vert") {
    return f_buff_vert;
  }
  if (name == "f_buff_frag") {
    return f_buff_frag;
  }

  if (name == "g_buff_vert") {
    return g_buff_vert;
  }
  if (name == "g_buff_frag") {
    return g_buff_frag;
  }

  if (name == "sqr_vert") {
    return sqr_vert;
  }
  if (name == "hdr_sqr_frag") {
    return hdr_sqr_frag;
  }

  if (name == "def_sqr_frag") {
    return def_sqr_frag;
  }

  if (name == "ssao_sqr_frag") {
    // find and replace SCREEN_WIDTH and SCREEN_HEIGHT
    std::string s = ssao_sqr_frag;
    std::string s1 = "SCREEN_WIDTH";
    std::string s2 = "SCREEN_HEIGHT";
    s.replace(s.find(s1), s1.length(), std::to_string(width));
    s.replace(s.find(s2), s2.length(), std::to_string(height));

    std::cout << name << std::endl;
    std::cout << s << std::endl;
    return s;
  }

  if (name == "bleed_sqr_frag") {
    std::string s = bleed_sqr_frag;
    std::string s1 = "SCREEN_WIDTH";
    std::string s2 = "SCREEN_HEIGHT";
    s.replace(s.find(s1), s1.length(), std::to_string(width));
    s.replace(s.find(s2), s2.length(), std::to_string(height));
    // dump name and source
    std::cout << name << std::endl;
    std::cout << s << std::endl;
    return s;
  }

  if (name == "blur_sqr_frag") {
    return blur_sqr_frag;
  }

  if (name == "dbg_point_vert") {
    return dbg_point_vert;
  }
  if (name == "dbg_point_frag") {
    return dbg_point_frag;
  }
  if (name == "dbg_point_geom") {
    return dbg_point_geo;
  }

  if (name == "dbg_line_vert") {
    return dbg_line_vert;
  }
  if (name == "dbg_line_frag") {
    return dbg_line_frag;
  }
  if (name == "dbg_line_geom") {
    return dbg_line_geo;
  }

  return "";
}
} // namespace gg
#endif