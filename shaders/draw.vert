#version 460

uniform float u_time;
uniform float u_draw_radius;

uniform mat4 u_projection;
uniform mat4 u_view;
uniform mat4 u_model;

layout(location = 0) in vec2 i_position;
layout(location = 1) in float i_phi;
layout(location = 2) in float i_n;
layout(location = 3) in float i_bin_index;

out VS_OUT
{
    vec3 color;
} vs_out;

vec3 hsv_to_rgb(in vec3 c) 
{
  const vec4 k = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  const vec3 p = abs(fract(c.xxx + k.xyz) * 6.0 - k.www);

  return c.z * mix(k.xxx, clamp(p - k.xxx, 0.0, 1.0), c.y);
}

float map(float value, float min1, float max1, float min2, float max2) 
{
  float final = min2 + (value - min1) * (max2 - min2) / (max1 - min1);
  return clamp(final, min2, max2);
}

const float pi = 3.14159;

void main() 
{
    gl_PointSize = max((pow(abs(i_n), 0.85)) * u_draw_radius, 2.0);
    //gl_PointSize = (1.0 - pow(abs(sin(pi + i_n * 0.5)), 0.5)) * u_draw_radius;

    gl_Position =  u_projection * u_view * u_model * vec4(i_position, 0.0, 1.0);

    vs_out.color = hsv_to_rgb(vec3(i_n, 1.0, 1.0));
    //vs_out.color = mix(vs_out.color, vec3(i_n, i_n * 0.33, 0.0), 0.6);

    vs_out.color = vec3(i_n, i_n * 0.33, 0.50);

    //vs_out.color = hsv2rgb(vec3(i_bin_index / 2500, 1.0, 1.0));
}