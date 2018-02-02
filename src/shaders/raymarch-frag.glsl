#version 300 es

precision highp float;

#define M_PI 3.1415926535897932384626433832795

out vec4 out_Col;

uniform mat4 u_ViewInv;
uniform vec2 u_Window;
uniform vec4 u_Time;

vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d )
{
    return a + b*cos( 6.28318*(c*t+d) );
}

float fract_comp(float x) {
    return x - floor(x);
}

float fade (float t) {
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); 
}

float lerp(float x0, float x1, float t) {
    return (1.0 - t) * x0 + (t * x1);
}

vec4 mod289(vec4 x)
{
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}
 
vec4 permute(vec4 x)
{
    return mod289(((x*34.0)+1.0)*x);
}

vec3 noise_gen3D(vec3 pos) {
    float x = fract(sin(dot(vec3(pos.x,pos.y,pos.z), vec3(12.9898, 78.233, 78.156))) * 43758.5453);
    float y = fract(sin(dot(vec3(pos.x,pos.y,pos.z), vec3(2.332, 14.5512, 170.112))) * 78458.1093);
    float z = fract(sin(dot(vec3(pos.x,pos.y,pos.z), vec3(400.12, 90.5467, 10.222))) * 90458.7764);
    return 2.0 * (vec3(x,y,z) - 0.5);
}

float dotGridGradient(vec3 grid, vec3 pos) {
    vec3 grad = normalize(noise_gen3D(grid));
    vec3 diff = (pos - grid);
    //return grad.x;
    return clamp(dot(grad,diff),-1.0,1.0);
}

float perlin3D(vec3 pos, float step) {
    pos = pos/step;
    vec3 ming = floor(pos / step) * step;
    ming = floor(pos);
    vec3 maxg = ming + vec3(step, step, step);
    maxg = ming + vec3(1.0);
    vec3 range = maxg - ming;
    vec3 diff = pos - ming;
    vec3 diff2 = maxg - pos;
    float d000 = dotGridGradient(ming, pos);
    float d001 = dotGridGradient(vec3(ming[0], ming[1], maxg[2]), pos);
    float d010 = dotGridGradient(vec3(ming[0], maxg[1], ming[2]), pos);
    float d011 = dotGridGradient(vec3(ming[0], maxg[1], maxg[2]), pos);
    float d111 = dotGridGradient(vec3(maxg[0], maxg[1], maxg[2]), pos);
    float d100 = dotGridGradient(vec3(maxg[0], ming[1], ming[2]), pos);
    float d101 = dotGridGradient(vec3(maxg[0], ming[1], maxg[2]), pos);
    float d110 = dotGridGradient(vec3(maxg[0], maxg[1], ming[2]), pos);

    float ix00 = mix(d000,d100, fade(diff[0]));
    float ix01 = mix(d001,d101, fade(diff[0]));
    float ix10 = mix(d010,d110, fade(diff[0]));
    float ix11 = mix(d011,d111, fade(diff[0]));

    float iy0 = mix(ix00, ix10, fade(diff[1]));
    float iy1 = mix(ix01, ix11, fade(diff[1]));

    float iz = mix(iy0, iy1, fade(diff[2]));

    //return abs(range.x * 9.0);
    return (iz + 1.0) / 2.0;
    //return range.y;
    //return abs(diff.x);
}

mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    angle = -angle;
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

vec3 opRt( vec3 p, vec3 axis, float angle )
{
    mat4 m = rotationMatrix(axis,angle);
    vec3 q = vec3(m*vec4(p,1.0));
    return q;
}

float opS( float d1, float d2 )
{
    return max(-d2,d1);
}
vec2 opSvec( vec2 d1, vec2 d2 )
{
    if(-d2.x > d1.x) {
        return vec2(-d2.x,d2.y);
    } else {
        return d1;
    }
}

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

float opI( float d1, float d2 )
{
    return max(d1,d2);
}

vec2 smin( vec2 a, vec2 b, float k )
{
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    float t = mix( b.x, a.x, h ) - k*h*(1.0-h);
    vec2 ret = vec2(t,a.y);
    return ret;
}

vec4 opUvec4(vec4 d1, vec4 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

vec4 sminvec4( vec4 a, vec4 b, float k )
{
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    float t = mix( b.x, a.x, h ) - k*h*(1.0-h);
    vec4 ret = vec4(t,a.yzw);
    return ret;
}


float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdEllipsoid( in vec3 p, in vec3 r )
{
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sdTorus( vec3 p, vec2 t )
{
    return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}

float sdCone( vec3 p, vec2 c )
{
    // c must be normalized
    p = opRt(p, vec3(1.0,0.0,0.0), -M_PI/2.0);
    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}

vec3 getRotatedPos(vec3 pos, vec3 rot) {
    float t = (cos(2.0 *M_PI * u_Time.y / (300.0 ) ) + 1.0 ) / 2.0;
    return opRt( pos, rot, M_PI/6.0 * (t - 0.5)  );
}

vec2 mapKoffing(vec3 pos, vec3 rot) {
    pos = getRotatedPos(pos,rot);
    vec2 sphere = vec2( opS(sdSphere(pos, 0.5),sdSphere(pos - vec3(0.0,0.0,1.0), 0.2) ), 0.0 );
    vec2 torus = vec2(sdTorus( opRt( pos, vec3(0.0,0.0,1.0), 0.0  ) - vec3(0.0,0.53,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 torus2 = vec2(sdTorus( opRt( pos, vec3(0.0,-0.2,1.0), M_PI/1.5  ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 torus3 = vec2(sdTorus( opRt( pos, vec3(1.0,-0.4,0.0), -M_PI/1.5  ) - vec3(0.0,0.5,0.0), vec2(0.04,0.015) ), 0.0 );
    vec2 torus4 = vec2(sdTorus( opRt( pos, vec3(0.0,0.4,1.0), M_PI/4.0 ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 torus5 = vec2(sdTorus( opRt( pos, vec3(0.0,-0.4,1.0), -M_PI/3.0 ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 torus6 = vec2(sdTorus( opRt( pos, vec3(0.0,-0.2,1.0), -M_PI/1.6  ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 torus7 = vec2(sdTorus( opRt( pos, vec3(1.0,0.0,0.0), -M_PI/6.0  ) - vec3(0.0,-0.53,0.0), vec2(0.07,0.02) ), 0.0 );
    vec2 body = smin(sphere,torus, 0.1);
    body = smin(body,torus2, 0.05);
    body = smin(body,torus3, 0.05);
    body = smin(body,torus4, 0.05);
    body = smin(body,torus5, 0.05);
    body = smin(body,torus6, 0.05);
    body = smin(body,torus7, 0.05);
    vec2 eye = vec2(opS(sdSphere(pos - vec3(-0.14,0.2,0.33), 0.13), sdBox(pos - vec3(-0.15,-0.28,0.33), vec3(0.5) )      ), 0.25 );
    eye = opU(eye, vec2(sdSphere(pos - vec3(-0.155,0.29,0.42), 0.01   ), 1.0 ));
    vec2 eye2 = vec2(opS(sdSphere(pos - vec3(0.14,0.2,0.33), 0.13), sdBox(pos - vec3(0.15,-0.28,0.33), vec3(0.5) )      ), 0.25);
    eye2 = opU(eye2, vec2(sdSphere(pos - vec3(0.155,0.29,0.42), 0.01   ), 1.0 ));
    eye = opU(eye2,eye);
    vec2 mouth = vec2(opS(sdSphere(pos - vec3(0.0,0.27,0.45), 0.35), sdBox(pos - vec3(0.0,1.12,0.5), vec3(1.0) )   ), 0.9);
    body = opSvec(body,mouth);
    vec2 lip = vec2(sdTorus(pos - vec3(0.0,0.12,0.17), vec2(0.3,0.02) ), 0.0 );
    body = smin(body,lip, 0.1);
    vec2 tooth = vec2(sdEllipsoid(pos - vec3(-0.22, 0.1,0.41), vec3(0.03,0.05,0.03) ), 0.25 );
    vec2 tooth2 = vec2(sdEllipsoid(pos - vec3(0.22, 0.1,0.41), vec3(0.03,0.05,0.03) ), 0.25 );
    body = opU(body,tooth);
    body = opU(body,tooth2);
    vec2 skull = vec2(sdTorus( opRt( pos * vec3(0.7,1.0,1.0), vec3(1.0,0.0,0.0), M_PI/3.0  ) - vec3(0.0,-0.46,0.0), vec2(0.06,0.04) ), 0.5 );
    vec2 capsule = vec2(sdCapsule(opRt( pos * vec3(0.7,1.0,1.0), vec3(1.0,0.0,0.0), M_PI/5.5  ), vec3(-0.1, -0.48, 0.1), vec3(0.1,-0.48,-0.1), 0.035 ), 0.5);
    vec2 capsule2 = vec2(sdCapsule(opRt( pos * vec3(0.7,1.0,1.0), vec3(1.0,0.0,0.0), M_PI/5.5  ), vec3(0.1, -0.48, 0.1), vec3(-0.1,-0.48,-0.1), 0.035 ), 0.5);
    capsule = smin(capsule,capsule2,0.01);
    skull = opU(skull,capsule);
    body = opU(body, skull);
    vec2 koffing = opU(body,eye);
    return koffing;
}

vec2 opRepKoffing( vec3 p, vec3 c )
{
    vec3 q = mod(p,c)-0.5*c;
    vec3 rot = floor(p/c);
    rot = vec3(noise_gen3D(rot));
    //q += vec3(0.0,0.0,5.0);
    
    return mapKoffing( q ,normalize(rot));
}

vec2 map(vec3 pos) {
	//return vec2(sdSphere(  pos - vec3(0.0,0.0,0.0)  , 1.0  ) , 1.0);
    
    vec2 res = mapKoffing(pos, vec3(1.0,1.0,0.0));
    if(u_Time.w > 0.0) {
        vec2 rep = opRepKoffing(pos, vec3(2.0,2.0,2.0));
        rep = vec2(opI(rep.x, sdBox(pos - vec3(0.0,0.0,-5.0),vec3(4.0,4.0,3.4)) ), rep.y);
        res = opU(res,rep);
    }

    return res;

}

vec4 map_smoke(vec3 pos) {
    pos = getRotatedPos(pos, vec3(1.0,1.0,0.0));
    //vec4 capsule1 = vec4(sdCapsule(pos, vec3(0.0,0.5,0.0), vec3(0.0,1.0,0.0), 0.05 ), vec3(0.0,1.0,0.0) );
    //vec4 capsule2 = vec4(sdCapsule(pos, vec3(0.5,0.5,0.0), vec3(1.0,1.0,0.0), 0.05 ), vec3(1.0,1.0,0.0) );
    float speed = 200.0;
    float up = fade( 1.0 * mod(u_Time.z,speed ) / speed);
    vec4 cone = vec4(opI(sdCone(pos - vec3(0.0,0.5,0.0), vec2(0.5,0.2)), sdBox(pos - vec3(0.0,up - 0.5,0.0), vec3(1.0)  )), vec3(0.0,0.5,0.0));
    vec4 ball = vec4(sdSphere(pos - vec3(0.0,up+0.5,0.0), 0.4 * up), vec3(0.0,0.5,0.0));
    cone = opUvec4(cone,ball);
    // vec3 rot1_pos = opRt( pos, vec3(1.0,0.0,0.0), M_PI/2.0  ) ;
    // vec3 rot1_dir = opRt( vec3(0.0,1.0,0.0), vec3(1.0,0.0,0.0), M_PI/2.0  ) ;
    // vec4 cone2 = vec4(opI(sdCone(rot1_pos - vec3(0.0,0.5,0.0), vec2(0.5,0.2)), sdBox(rot1_pos - vec3(0.0,0.5,0.0), vec3(1.0)  )), -normalize(rot1_dir) );

    float up2 = fade(1.0 * mod(u_Time.z + 23.0,speed ) / speed);
    vec3 rot1_pos = opRt( pos, vec3(0.0,-0.2,1.0), M_PI/1.5  ) ;
    vec3 rot1_dir = opRt( vec3(0.0,1.0,0.0), vec3(0.0,-0.2,1.0), -M_PI/1.5 ) ;
    vec4 cone2 = vec4(opI(sdCone(rot1_pos - vec3(0.0,0.5,0.0), vec2(0.5,0.2)), sdBox(rot1_pos - vec3(0.0,up2 - 0.5,0.0), vec3(1.0)  )), rot1_dir );
    vec4 ball2 = vec4(sdSphere(rot1_pos - vec3(0.0,up2+0.5,0.0), 0.4 * up2), rot1_dir);
    cone2 = opUvec4(cone2,ball2);



    float up3 = fade(1.0 * mod(u_Time.z + 82.0,speed ) / speed);

    vec2 torus6 = vec2(sdTorus( opRt( pos, vec3(0.0,-0.2,1.0), -M_PI/1.6  ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 0.0 );
    vec3 rot2_pos = opRt( pos, vec3(0.0,-0.2,1.0), -M_PI/1.6  ) ;
    vec3 rot2_dir = opRt( vec3(0.0,1.0,0.0), vec3(0.0,-0.2,1.0), M_PI/1.6 ) ;
    vec4 cone3 = vec4(opI(sdCone(rot2_pos - vec3(0.0,0.5,0.0), vec2(0.5,0.2)), sdBox(rot2_pos - vec3(0.0,up3 - 0.5,0.0), vec3(1.0)  )), rot2_dir );
    vec4 ball3 = vec4(sdSphere(rot2_pos - vec3(0.0,up3+0.5,0.0), 0.4 * up3), rot2_dir);
    cone3 = opUvec4(cone3,ball3);
    cone2 = opUvec4(cone2,cone3);
    return opUvec4(cone, cone2);
    //return opUvec4(capsule1,capsule2);
}

vec3 calcNormal(vec3 pos )
{
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy).x + 
					  e.yyx*map( pos + e.yyx).x + 
					  e.yxy*map( pos + e.yxy).x + 
					  e.xxx*map( pos + e.xxx).x );
}

//from iq
float calcAO(vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.06*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = map( aopos ).x;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}

vec3 calcSS(vec3 ray_Dir, vec3 vLight, vec3 norm, float thick) {
    vec3 vLTLight = vLight + norm;
    float fLTDot = pow(clamp(dot(-ray_Dir, -vLTLight), 0.0,1.0), 2.0) * 5.0;
    vec3 fLt = vec3(dot(norm, vLight) + dot(-ray_Dir, -vLight)) * (fLTDot + 0.2) * thick;
    return fLt;
}

float intersectSphere(vec3 rd, vec3 r0, vec3 s0, float sr){
    float a = dot(rd, rd);
    vec3 s0_r0 = r0 - s0;
    float b = 2.0 * dot(rd, s0_r0);
    float c = dot(s0_r0, s0_r0) - (sr * sr);
    if (b*b - 4.0*a*c < 0.0) {
        return -1.0;
    }
    return (-b - sqrt((b*b) - 4.0*a*c))/(2.0*a);
}

float fbm(vec3 P, vec3 scale) {
    float p1 = perlin3D(P, 0.5);
    float p2 = perlin3D(P, 0.2);
    float p3 = perlin3D(P, 0.1);
    float ret = p1 - (p2 * 0.7 ) + (p3 * 0.1);
    //ret *= 2.0;
    return clamp(ret,0.0,1.0);
}


vec3 velScale(vec3 vel) {
    vec3 stepdir = step(0.0,vel) - vec3(0.5);
    return vel * 2.0 + stepdir * 2.0;
}

vec2 raymarchSmoke(vec3 ray_O, vec3 ray_Dir, float t, float t_max, vec3 scale) {
    if(t_max <= 0.0) {t_max = t + 1.0;}
    int numSteps = 10;
    float step = 1.0/ float(numSteps);
    float t0 = t;
    float trans = 1.0;
    float lightEnergy = 0.0;

    vec3 lightDir = normalize(vec3(0,1,0));
	for (int i = 0; i < numSteps; i++) {
		vec3 P = ray_O + t * ray_Dir;
		//d += perlin3D(P, 0.1) * 0.05;
        float st = 0.0;
        //vec3 scale = vec3(2.0,2.0,1.0);
        //scale = normalize(P - scale);
        vec3 off_P = P - scale * 0.02 * u_Time.z;
        vec3 sP = P + st * lightDir;
        float dS = fbm(off_P + 0.1 * lightDir, scale );
        //if(length(P + 0.1 * lightDir) >= 1.0 ) dS = 0.0;
        float d = fbm(off_P, scale) * step;
        float absorbed = d;
        float shadow = smoothstep(1.0,0.0,dS) * 0.85 + 0.15;
        //shadow = 1.0;
        lightEnergy += absorbed * trans * shadow;
        trans *= 1.0 - d;
		t += step;
		if(t > t_max || trans < 0.01) {
            break;
		}
        // if(t > t_max || trans < 0.01 || length(P) > 1.001) {
        //     break;
		// }
	}
    lightEnergy *= 1.0;
	return vec2(lightEnergy,trans);
}

vec2 raymarchSdfs(vec3 ray_O, vec3 ray_Dir) {
	int count = 0;
	bool flag = true;
	float t = 0.0;

	while (flag) {
		vec3 P = ray_O + t * ray_Dir;
		vec2 map = map(P);
		float d = map.x;
		if(abs(d) < 0.001) {
			flag = false;
			return vec2(t, map.y);
		} 
		t += d;
		if(count > 126 || t > 100.0) {
			flag = false;
		}
		count++;
	}
	return vec2(0.0,-0.1);
}

vec4 raymarchSmokeBounds(vec3 ray_O, vec3 ray_Dir) {
	int count = 0;
	bool flag = true;
	float t = 0.0;

	while (flag) {
		vec3 P = ray_O + t * ray_Dir;
		vec4 map = map_smoke(P);
		float d = map.x;
		if(abs(d) < 0.001) {
			flag = false;
			return vec4(t, map.yzw);
		} 
		t += d;
		if(count > 126 || t > 200.0) {
			flag = false;
		}
		count++;
	}
	return vec4(-0.1);
}

vec3 getKoffingColor(float m) {
    vec3 col;
    if(m < 0.01) {
        col = vec3(0.6,0.6,0.8);
    } else if(m < 0.26) {
        col = vec3(1.0,1.0,1.0);
    } else if(m < 0.51) {
        col = vec3(0.92,0.91,0.8);
    } else if(m < 1.0) {
        col = vec3(0.77,0.54,0.61);
    } else {
        col = vec3(0.0);
    }
    return col;
}


void main() {
	float u = gl_FragCoord.x * 2.0 / u_Window.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / u_Window.y - 1.0;

    float imageWidth = u_Window.x;
    float imageHeight = u_Window.y;
    float imageAspectRatio = imageWidth / imageHeight; // assuming width > height 
    float fov = 45.0;
    float Px = (u) * tan(fov / 2.0 * M_PI / 180.0) * imageAspectRatio; 
    float Py = (v) * tan(fov / 2.0 * M_PI / 180.0);
    vec3 ray_O = vec3(0.0); 
    vec3 ray_Dir = vec3(Px, Py, -1) - ray_O;
    ray_Dir = normalize(ray_Dir);
    
    ray_O = vec3(u_ViewInv * vec4(ray_O , 1.0));
    ray_Dir = vec3(u_ViewInv * vec4(ray_Dir, 0.0));
	ray_Dir = normalize(ray_Dir);

	vec3 color = vec3(0.05,0.0,0.1);

	vec2 res = raymarchSdfs(ray_O, ray_Dir);
    vec3 lightvec = normalize(vec3(0.0,1.0,1.0));
    vec3 lightvecSS = normalize(vec3(0.0,0.0,-0.1));
	if(res.x > 0.0) {
		vec3 pos = ray_O + ray_Dir * res.x;
		vec3 norm = calcNormal(pos);
        float ao = calcAO(pos,norm);
        float thick = 1.0 -calcAO(pos,-norm);
        color = vec3(0.4,0.1,0.7);
        color = getKoffingColor(res.y);
        color *= clamp(ao + 0.2, 0.0,1.0);
        //color *= thick;
        
        vec3 sss = (calcSS(ray_Dir, lightvecSS, norm, thick)) ;
        color += color * sss * vec3(0.8,0.3,0.5);
        color *= clamp((dot(norm, lightvec) + 0.5)/1.5,0.0,1.0);
        color += clamp(dot(norm, lightvecSS) ,0.0,1.0) * vec3(0.8,0.2,0.4) * 0.2;
        vec3 H = normalize(lightvec + ray_Dir);
        color += pow(clamp(dot(H,norm),0.0,1.0 ),10.0) * vec3(1.0);
	}


    float smoke_bounds = intersectSphere(ray_Dir, ray_O, vec3(0,0,0), 1.0);
    smoke_bounds = -1.0;
	vec4 res_smoke = raymarchSmokeBounds(ray_O,ray_Dir);
    if (res_smoke.x > 0.0 ) {
        //color = res_smoke.yzw;
        smoke_bounds = res_smoke.x;
    }


    float smoke_bounds_masked = perlin3D(vec3(u * 5.0,v *5.0,ray_O.z + u_Time.z * 0.01), 0.8);
    smoke_bounds_masked *= clamp(perlin3D(vec3(u * 5.0,v *5.0,ray_O.z+ u_Time.z * 0.01), 0.2) * 0.8,0.0,1.0);
    //smoke_bounds_masked *= smoke_bounds;

    float v_density = 0.0;
    if(smoke_bounds > 0.0) {
        //color = vec3(1.0);
        vec3 smoke_scale = res_smoke.yzw;
        //color = normalize(smoke_scale);
        vec2 smoke_res = raymarchSmoke(ray_O,ray_Dir, smoke_bounds + 0.01, res.x , smoke_scale );
        v_density = smoke_res.x;
	    vec3 smoke_col = vec3(1.0, 1.0, 1.0);
        //color = mix(vec3(0.0,1.0,0.0), smoke_col, v_density);
        //smoke_col = mix(vec3(0.0,0.0,0.0), smoke_col, v_density);
        
        v_density = clamp((v_density - 0.1) * 2.0, 0.0, 1.0 );
        //smoke_col = pal(v_density, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
        smoke_col = mix(vec3(1.0,1.0,1.0) , vec3(0.8,0.8,0.0) , v_density );
        //v_density *= smoke_bounds_masked;
        if(smoke_bounds_masked < 0.1) {
            smoke_col = mix(color, smoke_col, smoke_bounds_masked * 10.0);
        }
        color = mix(color, smoke_col, (clamp((v_density),0.0,1.0) ) );
        
        // if(v_density < 0.1) {
        //     color = mix(vec3(0.05,0.0,0.1), color, v_density * 10.0 );
        // }
    }

    float screen_cloud = perlin3D(vec3(u * 5.0,v *5.0,ray_O.z), 0.9);

    screen_cloud *= perlin3D(vec3(u * 5.0,v *5.0,ray_O.z), 0.5);

    color = mix(vec3(0.6, 0.98, 0.6),color,(screen_cloud * 0.15) + 0.85 );
	out_Col = vec4(color, 1.0);
}
