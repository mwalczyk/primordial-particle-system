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

// https://observablehq.com/@daxreincarnate/primordial-particle-system#moveParticles

// Data that will be associated with the GLFW window
struct InputData
{
    bool imgui_active = false;
};

// Viewport and camera settings
const uint32_t window_w = 750;
const uint32_t window_h = 750;
bool first_mouse = true;
float last_x;
float last_y;
float zoom = 45.0f;
glm::mat4 arcball_camera_matrix = glm::lookAt(glm::vec3{ 0.0f, 0.0f, 4.0f }, glm::vec3{ 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.0f });
glm::mat4 arcball_model_matrix = glm::mat4{ 1.0f };

// Global settings
bool show_ui = true;
float draw_radius = 6.0;

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
    if (key == GLFW_KEY_U && action == GLFW_PRESS)
    {
        show_ui = !show_ui;
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

struct alignas(8) Particle
{
    glm::vec2 position;
    float phi;
    float n;
    float l;
    float r;
    float bin_index;

    glm::vec2 get_heading() const
    {
        return { cosf(phi), sinf(phi) };
    }
};

struct Params
{
    float alpha = 0.1f;
    float beta = 0.125f;
    float radius = 30.0f;
    float speed = 50.0f;
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

const size_t count = 5000;

class System
{

public:

    System(const Bounds& bounds, int32_t k = 4) :
        bounds{ bounds },
        k{ k }
    {
        reset();

        bin_size = 1 << k;
       // bin_size = 50;
        x_bins = static_cast<uint32_t>(ceilf((float)window_w / (float)bin_size));
        y_bins = static_cast<uint32_t>(ceilf((float)window_h / (float)bin_size));
        bins.resize(x_bins * y_bins);

        std::cout << "Bin size: " << bin_size << "\n";
        std::cout << "X bins: " << x_bins << "\n";
        std::cout << "Y bins: " << y_bins << "\n";

        initialize();
    }

    void set_params(const Params& updated)
    {
        params = updated;
    }

    void reset()
    {
        particles.clear();

        std::default_random_engine generator;
        std::uniform_real_distribution<float> distribution(0.0f, glm::pi<float>() * 2.0f);

        float max_dimension = std::max(window_w, window_h);

        for (size_t i = 0; i < count; ++i)
        {
            particles.push_back(Particle{
                glm::diskRand(max_dimension * 0.25f) + glm::vec2{window_w, window_h} * 0.5f,
                distribution(generator)
            });
        }
    }

    void initialize()
    {
        glCreateVertexArrays(1, &vao);

        glCreateBuffers(1, &vbo);
        glNamedBufferStorage(vbo, sizeof(Particle) * particles.size(), &particles[0], GL_DYNAMIC_STORAGE_BIT);

        {
            glVertexArrayVertexBuffer(vao, 0, vbo, 0, sizeof(Particle));

            // Enable vertex attributes
            glEnableVertexArrayAttrib(vao, 0);
            glEnableVertexArrayAttrib(vao, 1);
            glEnableVertexArrayAttrib(vao, 2);
            glEnableVertexArrayAttrib(vao, 3);

            // Set vertex attribute formats
            glVertexArrayAttribFormat(vao, 0, 2, GL_FLOAT, GL_FALSE, offsetof(Particle, position));
            glVertexArrayAttribFormat(vao, 1, 1, GL_FLOAT, GL_FALSE, offsetof(Particle, phi));
            glVertexArrayAttribFormat(vao, 2, 1, GL_FLOAT, GL_FALSE, offsetof(Particle, n));
            glVertexArrayAttribFormat(vao, 3, 1, GL_FLOAT, GL_FALSE, offsetof(Particle, bin_index));

            // Associate buffers with vertex attributes
            glVertexArrayAttribBinding(vao, 0, 0);
            glVertexArrayAttribBinding(vao, 1, 0);
            glVertexArrayAttribBinding(vao, 2, 0);
            glVertexArrayAttribBinding(vao, 3, 0);
        }
    }

    std::vector<Particle*> get_neighbors(const Particle& particle, float radius) 
    {
        std::vector<Particle*> region = get_region(
            static_cast<uint32_t>(particle.position.x - radius),
            static_cast<uint32_t>(particle.position.y - radius),
            static_cast<uint32_t>(particle.position.x + radius),
            static_cast<uint32_t>(particle.position.y + radius));

        std::vector<Particle*> neighbors;

        const float radius_squared = radius * radius;

        for (size_t i = 0; i < region.size(); ++i)
        {
            const Particle& current = *region[i];

            float delta_x = current.position.x - particle.position.x;
            float delta_y = current.position.y - particle.position.y;
            float distance_squared = delta_x * delta_x + delta_y * delta_y;

            if (distance_squared < radius_squared)
            {
                neighbors.push_back(region[i]);
            }
        }

        return neighbors;
    }

    std::vector<Particle*> get_region(uint32_t minX, uint32_t minY, uint32_t maxX, uint32_t maxY)
    {
        std::vector<Particle*> region;
        auto back = back_inserter(region);

        uint32_t min_x_bin = minX >> k;
        uint32_t max_x_bin = maxX >> k;
        uint32_t min_y_bin = minY >> k;
        uint32_t max_y_bin = maxY >> k;

        max_x_bin++;
        max_y_bin++;

        if (max_x_bin > x_bins)
        {
            max_x_bin = x_bins;
        }
        if (max_y_bin > y_bins)
        {
            max_y_bin = y_bins;
        }

        for (size_t y = min_y_bin; y < max_y_bin; ++y) 
        {
            for (size_t x = min_x_bin; x < max_x_bin; ++x)
            {
                std::vector<Particle*>& current = bins[y * x_bins + x];
                copy(current.begin(), current.end(), back);
            }
        }

        return region;
    }

    void update()
    {
        calculate_bins();
        step();
    }

    void calculate_bins()
    {
        for (auto& bin : bins)
        {
            bin.clear();
        }

        for (auto& particle : particles)
        {
            //auto bin_indices = glm::uvec2(particle.position) >> glm::uvec2{ k };

            const uint32_t x_bin = static_cast<uint32_t>(particle.position.x) >> k;
            const uint32_t y_bin = static_cast<uint32_t>(particle.position.y) >> k;
            const uint32_t bin_index = y_bin * x_bins + x_bin;

            particle.bin_index = bin_index;

            if (x_bin < x_bins && y_bin < y_bins)
            {
                bins[bin_index].push_back(&particle);
            }
        }
    }

    void step()
    {
        int32_t n_max = 0;

        for (auto& particle: particles)
        {   
            // Keep track of the number of left and right neighbors of `particle`
            int32_t l = 0;
            int32_t r = 0;

            // Get all neighbors of `particle`, i.e. the particles that reside in bins that are
            // less than `param.radius` units away
            auto neighbors = get_neighbors(particle, params.radius);

            for (const auto& other : neighbors)
            {
                if (other->position == particle.position)
                {
                    continue;
                }
                
                // Determine which side of `particle` this neighbor lies on using the cross product
                const auto a = other->position - particle.position;
                const auto b = particle.get_heading();
                const auto c = glm::cross(glm::vec3{ b, 0.0f }, glm::vec3{ a, 0.0f });

                if (c.z >= 0.0f)
                {
                    l++;
                }
                else
                {
                    r++;
                }
            }

            // Calculate the total number of neighbors, `n`
            const int32_t n = l + r;

            // Keep track of the particle that had the highest neighbor count
            if (n > n_max)
            {
                n_max = n;
            }

            particle.n = n;
            particle.l = l;
            particle.r = r;
        }
        
        for (auto& particle : particles)
        {
            // Change directions
            const float delta_phi = params.alpha + particle.n * params.beta * sign(particle.r - particle.l);
            particle.phi += delta_phi;

            // Calculate step direction + length
            auto velocity = glm::normalize(particle.get_heading());
            velocity *= params.speed * (particle.n / static_cast<float>(particles.size())); // 1.25f

            // Clamp to window bounds
            auto updated = particle.position + velocity;
            updated = glm::clamp(updated, glm::vec2{ 0.0f }, glm::vec2{ window_w, window_h });

            particle.position = updated;
        }

        for (auto& particle : particles)
        {
            particle.n /= n_max;
        }

        glNamedBufferSubData(vbo, 0, sizeof(Particle) * particles.size(), particles.data());
    }

    void draw() const
    {
        glBindVertexArray(vao);
        glDrawArrays(GL_POINTS, 0, particles.size());
    }

private:

    int32_t k;
    int32_t bin_size;
    int32_t x_bins;
    int32_t y_bins;
    std::vector<std::vector<Particle*>> bins;

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
    auto projection = glm::ortho(0.0f, static_cast<float>(window_w), 0.0f, static_cast<float>(window_h));

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            ImGui::Begin("Primordial Particle System");
            ImGui::SliderFloat("Alpha", &params.alpha, -0.5f, 0.5f);
            ImGui::SliderFloat("Beta", &params.beta, -0.5f, 0.5f);
            ImGui::SliderFloat("Radius", &params.radius, 5.0f, 50.0f);
            ImGui::SliderFloat("Speed", &params.speed, 10.0f, 100.0f);
            ImGui::SliderFloat("Draw Radius", &draw_radius, 1.0f, 20.0f);
            if (ImGui::Button("Reset"))
            {
                system.reset();
            }
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
            

            system.set_params(params);
            system.update();

            shader_draw.use();
            shader_draw.uniform_float("u_time", glfwGetTime());
            shader_draw.uniform_float("u_draw_radius", draw_radius);
            shader_draw.uniform_mat4("u_model", glm::mat4{ 1.0f });
            shader_draw.uniform_mat4("u_view", glm::mat4{ 1.0f });
            shader_draw.uniform_mat4("u_projection", projection);

            system.draw();
        }

        // Draw the ImGui window
        {
            if (show_ui)
            {
                ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            }
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

