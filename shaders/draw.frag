#version 460

layout(location = 0) out vec4 o_color;

in VS_OUT
{
    vec3 color;
} fs_in;

void main() 
{	
	vec2 coord = gl_PointCoord - vec2(0.5);  
	if (length(coord) > 0.5) discard;
	
	o_color = vec4(fs_in.color, 1.0);
}