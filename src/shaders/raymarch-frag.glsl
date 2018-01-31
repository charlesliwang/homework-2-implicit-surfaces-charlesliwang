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

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}



vec2 smin( vec2 a, vec2 b, float k )
{
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    float t = mix( b.x, a.x, h ) - k*h*(1.0-h);
    float m = mix( b.y, a.y, h );
    vec2 ret = vec2(t,m);
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

vec2 mapKoffing(vec3 pos) {
    vec2 sphere = vec2( opS(sdSphere(pos, 0.5),sdSphere(pos - vec3(0.0,0.0,1.0), 0.2) ), 1.0 );
    vec2 torus = vec2(sdTorus( opRt( pos, vec3(0.0,0.0,1.3), M_PI/6.0  ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 1.0 );
    vec2 torus2 = vec2(sdTorus( opRt( pos, vec3(0.0,0.0,1.3), -M_PI/1.5  ) - vec3(0.0,0.5,0.0), vec2(0.07,0.02) ), 1.0 );
    vec2 body = smin(sphere,torus, 0.1);
    body = smin(body,torus2, 0.1);
    vec2 eye = vec2(opS(sdSphere(pos - vec3(-0.14,0.2,0.33), 0.13), sdBox(pos - vec3(-0.15,-0.28,0.33), vec3(0.5) )      ));
    eye = opU(eye, vec2(sdSphere(pos - vec3(-0.155,0.29,0.42), 0.01   ), 1.0 ));
    vec2 eye2 = vec2(opS(sdSphere(pos - vec3(0.14,0.2,0.33), 0.13), sdBox(pos - vec3(0.15,-0.28,0.33), vec3(0.5) )      ));
    eye2 = opU(eye2, vec2(sdSphere(pos - vec3(0.155,0.29,0.42), 0.01   ), 1.0 ));
    // vec2 eye2 = vec2(opS(sdSphere(pos - vec3(0.14,0.13,0.36), 0.13), sdBox(pos - vec3(0.15,-0.37,0.35), vec3(0.5) )      ));
    // eye2 = opU(eye2, vec2(sdSphere(pos - vec3(0.155,0.22,0.45), 0.01   ), 1.0 ));
    eye = opU(eye2,eye);
    vec2 mouth = vec2(opS(sdSphere(pos - vec3(0.0,0.27,0.45), 0.35), sdBox(pos - vec3(0.0,1.12,0.5), vec3(1.0) )      ));
    body = vec2(opS(body.x,mouth.x),body.x );
    vec2 lip = vec2(sdTorus(pos - vec3(0.0,0.12,0.17), vec2(0.3,0.02) ), 1.0 );
    body = smin(body,lip, 0.1);
    vec2 tooth = vec2(sdEllipsoid(pos - vec3(-0.22, 0.1,0.41), vec3(0.03,0.05,0.03) ), 1.0 );
    vec2 tooth2 = vec2(sdEllipsoid(pos - vec3(0.22, 0.1,0.41), vec3(0.03,0.05,0.03) ), 1.0 );
    body = opU(body,tooth);
    body = opU(body,tooth2);
    vec2 koffing = opU(body,eye);
    return koffing;
}

vec2 opRepKoffing( vec3 p, vec3 c )
{
    vec3 q = mod(p,c)-0.5*c;
    return mapKoffing( q );
}

vec2 map(vec3 pos) {
	//return vec2(sdSphere(  pos - vec3(0.0,0.0,0.0)  , 1.0  ) , 1.0);
    //vec2 rep = opRepKoffing(pos, vec3(100.0));
    vec2 res = mapKoffing(pos);
    return res;

}

vec3 calcNormal(vec3 pos )
{
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy).x + 
					  e.yyx*map( pos + e.yyx).x + 
					  e.yxy*map( pos + e.yxy).x + 
					  e.xxx*map( pos + e.xxx).x );
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
    float p1 = perlin3D(P / scale, 0.8);
    float p2 = perlin3D(P / scale, 0.3);
    float p3 = perlin3D(P / scale, 0.1);
    float ret = p1 - (p2 * 1.0 ) + (p3 * 0.2);
    ret *= 2.0;
    return clamp(0.0,1.0,ret);
}

vec2 raymarchSmoke(vec3 ray_O, vec3 ray_Dir, float t, float t_max) {
    if(t_max <= 0.0) {t_max = t + 1.0;}
    int numSteps = 20;
    float step = 1.0/ float(numSteps);
    float t0 = t;
    float trans = 1.0;
    float lightEnergy = 0.0;

    vec3 lightDir = normalize(vec3(0,1,0));
	for (int i = 0; i < numSteps; i++) {
		vec3 P = ray_O + t * ray_Dir;
		//d += perlin3D(P, 0.1) * 0.05;
        float st = 0.0;
        vec3 scale = vec3(2.0,2.0,1.0);
        vec3 off_P = P - normalize(scale - vec3(1.0)) * 0.02 * u_Time.z;
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
		if(t > t_max || trans < 0.01 || length(P) > 1.001) {
            break;
		}
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
		if(count > 126 || t > 200.0) {
			flag = false;
		}
		count++;
	}
	return vec2(0.0,-0.1);
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
	if(res.x > 0.0) {
		vec3 pos = ray_O + ray_Dir * res.x;
		vec3 norm = calcNormal(pos);
        color = vec3(0.4,0.1,0.7);
        color *= clamp(0.0,1.0,dot(norm,vec3(0.0,1.0,0.0) ) );
        color = abs(norm);
	}
    float t_sphere = intersectSphere(ray_Dir, ray_O, vec3(0,0,0), 1.0);
    t_sphere = -1.0;
    float v_density = 0.0;
    if(t_sphere > 0.0) {
        //color = vec3(1.0);
        vec2 smoke_res = raymarchSmoke(ray_O,ray_Dir, t_sphere + 0.01, res.x );
        v_density = smoke_res.x;
	    vec3 smoke_col = vec3(1.0, 1.0, 1.0);
        //color = mix(vec3(0.0,1.0,0.0), smoke_col, v_density);
        //smoke_col = mix(vec3(0.0,0.0,0.0), smoke_col, v_density);
        
        v_density = clamp( 0.0, 1.0, (v_density - 0.2) / 0.7 );
        //smoke_col = pal(v_density, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
        
        smoke_col = mix(vec3(1.0,1.0,1.0) , vec3(0.8,0.8,0.0) , v_density + 0.5 );
        color = mix(color, smoke_col, (clamp(0.0,1.0,(v_density + 0.2)) ) );
        // if(v_density < 0.1) {
        //     color = mix(vec3(0.05,0.0,0.1), color, v_density * 10.0 );
        // }
    }

	out_Col = vec4(color, 1.0);
}
