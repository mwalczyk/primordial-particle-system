#include <algorithm>
#include <iostream>
#include <random>

#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/random.hpp"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "shader.h"

// Data that will be associated with the GLFW window
struct InputData
{
    bool imgui_active = false;
};

// Viewport and camera settings
const uint32_t window_w = 800;
const uint32_t window_h = 800;
bool first_mouse = true;
float last_x;
float last_y;
float zoom = 45.0f;
glm::mat4 arcball_camera_matrix = glm::lookAt(glm::vec3{ 0.0f, 0.0f, 4.0f }, glm::vec3{ 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.0f });
glm::mat4 arcball_model_matrix = glm::mat4{ 1.0f };

// Global settings
// ...

// Appearance settings
ImVec4 clear_color = ImVec4(0.66f, 0.64f, 0.66f, 1.0f);

InputData input_data;

/**
 * A function for handling scrolling.
 */
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (zoom >= 1.0f && zoom <= 45.0f)
    {
        zoom -= yoffset;
    }
    if (zoom <= 1.0f)
    {
        zoom = 1.0f;
    }
    if (zoom >= 45.0f)
    {
        zoom = 45.0f;
    }
}

/**
 * A function for handling key presses.
 */
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        // Close the GLFW window
        glfwSetWindowShouldClose(window, true);
    }
    if (key == GLFW_KEY_H && action == GLFW_PRESS)
    {
        // Reset the arcball camera
        arcball_camera_matrix = glm::lookAt(glm::vec3{ 6.0f, 0.0f, 0.0f }, glm::vec3{ 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.0f });
        arcball_model_matrix = glm::mat4{ 1.0f };
    }
}

/**
 * Get a normalized vector from the center of a virtual sphere centered at the origin to
 * a point `point_on_sphere` on the virtual ball surface, such that `point_on_sphere`
 * is aligned on screen's (x, y) coordinates.  If (x, y) is too far away from the
 * sphere, return the nearest point on the virtual ball surface.
 */
glm::vec3 get_arcball_vector(int x, int y)
{
    auto point_on_sphere = glm::vec3{
        1.0f * x / window_w * 2.0f - 1.0f,
        1.0f * y / window_h * 2.0f - 1.0f,
        0.0f
    };

    point_on_sphere.y = -point_on_sphere.y;

    const float op_squared = point_on_sphere.x * point_on_sphere.x + point_on_sphere.y * point_on_sphere.y;

    if (op_squared <= 1.0f * 1.0f)
    {
        // Pythagorean theorem
        point_on_sphere.z = sqrt(1.0f * 1.0f - op_squared);
    }
    else
    {
        // Nearest point
        point_on_sphere = glm::normalize(point_on_sphere);
    }

    return point_on_sphere;
}

/**
 * Performs arcball camera calculations.
 */
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    // First, check if the user is interacting with the ImGui interface - if they are,
    // we don't want to process mouse events any further
    auto input_data = static_cast<InputData*>(glfwGetWindowUserPointer(window));

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && !input_data->imgui_active)
    {
        if (first_mouse)
        {
            last_x = xpos;
            last_y = ypos;
            first_mouse = false;
        }

        if (xpos != last_x || ypos != last_y)
        {
            const float rotation_speed = 0.25f;

            glm::vec3 va = get_arcball_vector(last_x, last_y);
            glm::vec3 vb = get_arcball_vector(xpos, ypos);
            const float angle = acos(std::min(1.0f, glm::dot(va, vb))) * rotation_speed;
            const glm::vec3 axis_camera_coordinates = glm::cross(va, vb);

            glm::mat3 camera_to_object = glm::inverse(glm::mat3(arcball_camera_matrix) * glm::mat3(arcball_model_matrix));

            glm::vec3 axis_in_object_coord = camera_to_object * axis_camera_coordinates;

            arcball_model_matrix = glm::rotate(arcball_model_matrix, glm::degrees(angle), axis_in_object_coord);

            // Set last to current
            last_x = xpos;
            last_y = ypos;
        }
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE)
    {
        last_x = xpos;
        last_y = ypos;
    }
}

/**
 * Debug function that will be used internally by OpenGL to print out warnings, errors, etc.
 */
void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const* message, void const* user_param)
{
    auto const src_str = [source]() {
        switch (source)
        {
        case GL_DEBUG_SOURCE_API: return "API";
        case GL_DEBUG_SOURCE_WINDOW_SYSTEM: return "WINDOW SYSTEM";
        case GL_DEBUG_SOURCE_SHADER_COMPILER: return "SHADER COMPILER";
        case GL_DEBUG_SOURCE_THIRD_PARTY: return "THIRD PARTY";
        case GL_DEBUG_SOURCE_APPLICATION: return "APPLICATION";
        case GL_DEBUG_SOURCE_OTHER: return "OTHER";
        }
    }();

    auto const type_str = [type]() {
        switch (type)
        {
        case GL_DEBUG_TYPE_ERROR: return "ERROR";
        case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATED_BEHAVIOR";
        case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: return "UNDEFINED_BEHAVIOR";
        case GL_DEBUG_TYPE_PORTABILITY: return "PORTABILITY";
        case GL_DEBUG_TYPE_PERFORMANCE: return "PERFORMANCE";
        case GL_DEBUG_TYPE_MARKER: return "MARKER";
        case GL_DEBUG_TYPE_OTHER: return "OTHER";
        }
    }();

    auto const severity_str = [severity]() {
        switch (severity) {
        case GL_DEBUG_SEVERITY_NOTIFICATION: return "NOTIFICATION";
        case GL_DEBUG_SEVERITY_LOW: return "LOW";
        case GL_DEBUG_SEVERITY_MEDIUM: return "MEDIUM";
        case GL_DEBUG_SEVERITY_HIGH: return "HIGH";
        }
    }();

    std::cout << src_str << ", " << type_str << ", " << severity_str << ", " << id << ": " << message << '\n';
}

struct Bounds
{
    glm::vec2 min;
    glm::vec2 max;
};

struct Particle
{
    glm::vec2 position;
    float phi;
    float n;

    glm::vec2 get_heading() const
    {
        return { cosf(phi), sinf(phi) };
    }
};

struct Params
{
    float alpha = 0.1f;
    float beta = 0.5f;
    float radius = 0.05f;
    float speed = 0.0015f;
};

/*
 * Returns the sign of `val`. 
 *
 * Referenced from: https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
 */
template <typename T> int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}

class System
{

public:

    System(const Bounds& bounds, size_t number_of_particles = 2000) :
        bounds{ bounds }
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<float> distribution(0.0f, glm::pi<float>() * 2.0f);

        for (size_t i = 0; i < number_of_particles; ++i)
        {
            particles.push_back(Particle{ glm::diskRand(0.5f) + 0.5f, distribution(generator) });
        }

        initialize();
    }

    void set_params(const Params& updated)
    {
        params = updated;
    }

    void initialize()
    {
        glCreateVertexArrays(1, &vao);

        glCreateBuffers(1, &vbo);
        glNamedBufferStorage(vbo, sizeof(Particle) * particles.size(), &particles[0], GL_DYNAMIC_STORAGE_BIT);

        {
            glVertexArrayVertexBuffer(vao, 0, vbo, 0, sizeof(Particle));

            glEnableVertexArrayAttrib(vao, 0);
            glEnableVertexArrayAttrib(vao, 1);
            glEnableVertexArrayAttrib(vao, 2);

            glVertexArrayAttribFormat(vao, 0, 2, GL_FLOAT, GL_FALSE, offsetof(Particle, position));
            glVertexArrayAttribFormat(vao, 1, 1, GL_FLOAT, GL_FALSE, offsetof(Particle, phi));
            glVertexArrayAttribFormat(vao, 2, 1, GL_FLOAT, GL_FALSE, offsetof(Particle, n));

            glVertexArrayAttribBinding(vao, 0, 0);
            glVertexArrayAttribBinding(vao, 1, 0);
            glVertexArrayAttribBinding(vao, 2, 0);
        }
    }

    void update()
    {
        int32_t n_max = 0;

        for (size_t i = 0; i < particles.size(); ++i)
        {
            int32_t l = 0;
            int32_t r = 0;

            for (size_t j = 0; j < particles.size(); ++j)
            {
                if (i != j)
                {
                    if (glm::distance(particles[j].position, particles[i].position) < params.radius)
                    {
                        auto a = particles[j].position - particles[i].position;
                        auto b = particles[i].get_heading();

                        auto c = glm::cross(glm::vec3{ b, 0.0f }, glm::vec3{ a, 0.0f });

                        if (c.z >= 0.0f)
                        {
                            l++;
                        }
                        else
                        {
                            r++;
                        }
                    }

                }
            }

            int32_t n = l + r;

            if (n > n_max)
            {
                n_max = n;
            }

            // Change directions
            const float delta_phi = params.alpha + n * params.beta * sign(r - l);
            particles[i].phi += delta_phi;

            // Calculate step direction + length
            const float max_speed = 0.25f;
            auto velocity = glm::normalize(particles[i].get_heading());
            velocity *= max_speed * (n / static_cast<float>(particles.size()));// params.speed;
            
            // Clamp to window bounds
            auto updated = particles[i].position + velocity;
            updated = glm::clamp(updated, glm::vec2{ 0.0f }, glm::vec2{ 1.0f });

            particles[i].position = updated;
            particles[i].n = n;
        }
        
        for (auto& particle : particles)
        {
            particle.n /= n_max;
        }

        glNamedBufferSubData(vbo, 0, sizeof(Particle) * particles.size(), &particles[0]);
    }

    void draw() const
    {
        glBindVertexArray(vao);
        glDrawArrays(GL_POINTS, 0, particles.size());
    }

private:

    Bounds bounds;
    std::vector<Particle> particles;
    Params params;
    uint32_t vao;
    uint32_t vbo;

};

int main()
{
    // Create and configure the GLFW window 
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, false);
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow* window = glfwCreateWindow(window_w, window_h, "Marching Cubes", nullptr, nullptr);

    if (window == nullptr)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetWindowUserPointer(window, &input_data);

    // Load function pointers from glad
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");

    // Setup initial OpenGL state
    {
#if defined(_DEBUG)
        // Debug logging
        glEnable(GL_DEBUG_OUTPUT);
        glDebugMessageCallback(message_callback, nullptr);
#endif
        // Let the vertex shader set the point size
        glEnable(GL_PROGRAM_POINT_SIZE);
    }

    auto bounds = Bounds{ glm::vec2{0.0f}, glm::vec2{1.0f} };
    auto params = Params{};
    auto system = System{ bounds };

    auto shader_draw = graphics::Shader{ "../shaders/draw.vert", "../shaders/draw.frag" };

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            ImGui::Begin("Primordial Particle System");
            ImGui::SliderFloat("Alpha", &params.alpha, -glm::pi<float>() * 0.25f, glm::pi<float>() * 0.25f);
            ImGui::SliderFloat("Beta", &params.beta, -glm::pi<float>() * 0.25f, glm::pi<float>() * 0.25f);
            ImGui::SliderFloat("Radius", &params.radius, 0.001f, 0.2f);
            ImGui::SliderFloat("Speed", &params.speed, 0.001f, 0.01f);
            ImGui::ColorEdit3("clear color", (float*)&clear_color);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }
        ImGui::Render();

        glViewport(0, 0, window_w, window_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw the particle system
        {
            auto projection = glm::ortho(-0.5f, 1.5f, -0.5f, 1.5f);

            system.set_params(params);
            system.update();

            shader_draw.use();
            shader_draw.uniform_float("u_time", glfwGetTime());
            shader_draw.uniform_mat4("u_model", glm::mat4{ 1.0f });
            shader_draw.uniform_mat4("u_view", glm::mat4{ 1.0f });
            shader_draw.uniform_mat4("u_projection", projection);

            system.draw();
        }

        // Draw the ImGui window
        {
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        glfwSwapBuffers(window);
    }

    // Clean-up imgui
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
       
    // Clean-up GLFW
    glfwDestroyWindow(window);
    glfwTerminate();
}

