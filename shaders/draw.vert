#version 460

uniform float u_time;

uniform mat4 u_projection;
uniform mat4 u_view;
uniform mat4 u_model;

layout(location = 0) in vec2 i_position;
layout(location = 1) in float i_phi;
layout(location = 2) in float i_n;

out VS_OUT
{
    vec3 color;
} vs_out;

vec3 hsv2rgb(vec3 c) 
{
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

float map(float value, float min1, float max1, float min2, float max2) 
{
  return min2 + (value - min1) * (max2 - min2) / (max1 - min1);
}

const float pi = 3.14159;

void main() 
{
   // gl_PointSize = max(pow(i_n, 0.525) * 6.0, 2.0);
    gl_PointSize = max(pow(i_n, 2.0) * 10.0, 2.0);
    gl_PointSize = (pow(abs(i_n), 0.85)) * 8.0;
    //gl_PointSize = (1.0 - pow(abs(sin(pi + i_n * 0.5)), 0.5)) * 12.0;

    gl_Position =  u_projection * u_view * u_model * vec4(i_position, 0.0, 1.0);

    vs_out.color = hsv2rgb(vec3(map(i_n, 0.0, 1.0, 0.15, 0.85), 1.0, 1.0));
    vs_out.color = mix(vs_out.color, vec3(i_n, i_n * 0.33, 0.0), 0.5);
}