#include "deformable_mesh.h"

#include <GLFW/glfw3.h>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <igl/file_dialog_open.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <random>

struct pick_t
{
    unsigned int k  = 0;
    unsigned int vi = 0;
    bool picked{false};
    mti855::physics::deformable_mesh* object;
    Eigen::Vector3f bc;
};

struct simulation_params_t
{
    float Rb              = 0.1f;
    float beta            = 0.8f;
    float Famplitude      = 10.f;
    float PickAmplifier   = 0.1f;
    float dt              = 0.01f;
    float tau             = 0.8f;
    float perturbation    = 0.1f;
    bool visualize_forces = false;
    bool show_triangles   = false;
    mti855::physics::deformation_type_t deformation_type =
        mti855::physics::deformation_type_t::linear;
    bool draw_bounding_box        = false;
    bool has_gravity              = false;
    bool pause                    = false;
    unsigned int fps              = 0;
    unsigned int vertices         = 0;
    unsigned int triangles        = 0;
    unsigned int total_memory     = 0;
    unsigned int bytes_per_vertex = 0;
};

int main(int argc, char* argv[])
{
    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<> rd(0.f, 1.f);

    igl::opengl::glfw::Viewer viewer;

    mti855::physics::deformable_mesh object;
    std::atomic<bool> object_loaded{false};

    viewer.data().show_labels = true;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    simulation_params_t sim_params;

    menu.callback_draw_viewer_window = [&]() {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(10.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 600), ImGuiCond_FirstUseEver);
        ImGui::Begin("Meshless Deformations", nullptr, ImGuiWindowFlags_NoSavedSettings);

        bool dirty = false;

        if (ImGui::CollapsingHeader("Model", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Load OBJ Triangle Mesh"))
            {
                std::string const filename = igl::file_dialog_open();
                std::filesystem::path const mesh{filename};
                if (std::filesystem::exists(mesh) && std::filesystem::is_regular_file(mesh))
                {
                    object.load_from_obj(filename);
                    object_loaded = true;
                    dirty         = true;
                }
            }
        }

        // Expose the same variable directly ...
        if (ImGui::CollapsingHeader("Integration", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::SliderFloat("Velocity Damping", &sim_params.Rb, 0.f, 10.f);
            ImGui::SliderFloat("Force Amplitude", &sim_params.Famplitude, 0.f, 1000.f);
            ImGui::SliderFloat("Picking Amplifier", &sim_params.PickAmplifier, 0.f, 1.f);
            ImGui::SliderFloat("Tau", &sim_params.tau, 0.f, 1.f);
            ImGui::SliderFloat("Beta", &sim_params.beta, 0.f, 1.f);
            ImGui::SliderFloat("Regularization Perturbation", &sim_params.perturbation, 0.f, 0.1f);
            ImGui::SliderFloat("Time step", &sim_params.dt, 0.001f, 1.0f);
            ImGui::Checkbox("Activate gravity", &sim_params.has_gravity);
            ImGui::Checkbox("Pause", &sim_params.pause);
        }

        if (ImGui::CollapsingHeader("Deformation types", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::RadioButton("Rotation", (int*)(&sim_params.deformation_type), 0))
            {
                sim_params.deformation_type = mti855::physics::deformation_type_t::rotation;
            }
            if (ImGui::RadioButton("Linear", (int*)(&sim_params.deformation_type), 1))
            {
                sim_params.deformation_type = mti855::physics::deformation_type_t::linear;
            }
            if (ImGui::RadioButton("Quadratic", (int*)(&sim_params.deformation_type), 2))
            {
                sim_params.deformation_type = mti855::physics::deformation_type_t::quadratic;
            }
        }

        if (ImGui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Show force field", &sim_params.visualize_forces);
            ImGui::Checkbox("Show triangles", &sim_params.show_triangles);
            ImGui::Checkbox("Show bounding box", &sim_params.draw_bounding_box);
        }

        if (ImGui::CollapsingHeader("Data", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Text("FPS: %d", sim_params.fps);
            ImGui::Text("Vertices: %d", sim_params.vertices);
            ImGui::Text("Triangles: %d", sim_params.triangles);
            ImGui::Text("Total Memory (Kb): %d", sim_params.total_memory / 1024);
            ImGui::Text(
                "Bytes per vertex: %d",
                sim_params.vertices > 0u ? sim_params.total_memory / sim_params.vertices : 0u);
        }

        if (dirty)
        {
            viewer.data().clear();
            viewer.data().set_mesh(object.V(), object.F());
            viewer.core().align_camera_center(object.V());
        }

        ImGui::End();
    };

    viewer.core().is_animating = true;

    std::atomic<bool> new_physics_update = false;

    pick_t pick;
    viewer.callback_mouse_down =
        [&pick, &object](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        using button_type = igl::opengl::glfw::Viewer::MouseButton;
        if (static_cast<button_type>(button) != button_type::Left)
            return false;

        if (modifier == GLFW_MOD_SHIFT)
        {
            double const x = static_cast<double>(viewer.current_mouse_x);
            double const y =
                viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);

            int fid;
            Eigen::Vector3f bc{};

            if (igl::unproject_onto_mesh(
                    Eigen::Vector2f(x, y),
                    viewer.core().view,
                    viewer.core().proj,
                    viewer.core().viewport,
                    object.V(),
                    object.F(),
                    fid,
                    bc))
            {
                auto const& F = object.F();
                Eigen::Vector3i const face{F(fid, 0), F(fid, 1), F(fid, 2)};
                unsigned int closest_vertex = face(0);

                if (bc(1) > bc(0) && bc(1) > bc(2))
                {
                    closest_vertex = face(1);
                }
                else if (bc(2) > bc(0) && bc(2) > bc(1))
                {
                    closest_vertex = face(2);
                }

                object.set_fixed(closest_vertex, !object.is_fixed(closest_vertex));
            }

            return true;
        }

        if (modifier == GLFW_MOD_CONTROL)
        {
            // do picking
            int fid;
            double const x = static_cast<double>(viewer.current_mouse_x);
            double const y =
                viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);

            Eigen::Vector3f bc{};

            bool const hit = igl::unproject_onto_mesh(
                Eigen::Vector2f(x, y),
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport,
                object.V(),
                object.F(),
                fid,
                bc);

            if (hit)
            {
                pick.picked = true;
                pick.k      = fid;
                pick.object = &object;
                pick.bc     = bc;

                auto const& F = pick.object->F();
                Eigen::Vector3i const face{F(pick.k, 0), F(pick.k, 1), F(pick.k, 2)};
                unsigned int closest_vertex = face(0);

                if (pick.bc(1) > pick.bc(0) && pick.bc(1) > pick.bc(2))
                {
                    closest_vertex = face(1);
                }
                else if (pick.bc(2) > pick.bc(0) && pick.bc(2) > pick.bc(1))
                {
                    closest_vertex = face(2);
                }

                pick.vi = closest_vertex;
            }

            return true;
        }

        return false;
    };

    viewer.callback_mouse_up =
        [&pick](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        using button_type = igl::opengl::glfw::Viewer::MouseButton;
        if (static_cast<button_type>(button) != button_type::Left)
            return false;

        if (pick.picked)
        {
            pick.picked = false;
        }

        return true;
    };

    enum class user_input_force_t { none, up, down, left, right, forward, backward };

    std::atomic<user_input_force_t> user_input_force = user_input_force_t::none;
    viewer.callback_key_pressed                      = [&user_input_force](
                                      igl::opengl::glfw::Viewer& viewer,
                                      unsigned int key,
                                      int modifiers) -> bool {
        switch (key)
        {
            case 'i': user_input_force = user_input_force_t::up; break;
            case 'k': user_input_force = user_input_force_t::down; break;
            case 'j': user_input_force = user_input_force_t::left; break;
            case 'l': user_input_force = user_input_force_t::right; break;
            case 'b': user_input_force = user_input_force_t::backward; break;
            case GLFW_KEY_SPACE: user_input_force = user_input_force_t::forward; break;
            default: return false;
        }

        return true;
    };

    auto const& simulate = [&]() {
        if (!object_loaded)
            return;

        if (user_input_force != user_input_force_t::none)
        {
            Eigen::Vector3d const geometric_center = object.V().colwise().mean();
            Eigen::Vector4d force4d{0., 0., 0., 1.};

            switch (user_input_force)
            {
                case user_input_force_t::up: force4d(1) = 1.; break;
                case user_input_force_t::down: force4d(1) = -1.; break;
                case user_input_force_t::left: force4d(0) = -1.; break;
                case user_input_force_t::right: force4d(0) = 1.; break;
                case user_input_force_t::forward: force4d(2) = 1.; break;
                case user_input_force_t::backward: force4d(2) = -1.; break;
            }

            force4d *= sim_params.Famplitude;
            force4d.w() = 1.0;

            // convert force to world space
            Eigen::Matrix4d const screen_to_world_transform =
                (viewer.core().proj * viewer.core().view).inverse().cast<double>();
            Eigen::Vector3d const force = (screen_to_world_transform * force4d).segment(0, 3);

            Eigen::MatrixXd const f = object.apply_force(geometric_center, force);

            if (sim_params.visualize_forces)
            {
                Eigen::VectorXd const s = f.rowwise().norm();
                viewer.data().set_data(s);
            }

            user_input_force = user_input_force_t::none;
        }

        if (sim_params.has_gravity)
        {
            object.apply_gravity();
        }

        if (!sim_params.visualize_forces)
        {
            Eigen::VectorXd s(object.V().rows());
            s.setConstant(1.0);
            viewer.data().set_data(s);
        }

        object.set_deformation_type(sim_params.deformation_type);
        object.set_tau(sim_params.tau);
        object.set_rayleigh_beta(sim_params.Rb);
        object.set_beta(sim_params.beta);
        object.set_perturbation(sim_params.perturbation);

        if (sim_params.deformation_type == mti855::physics::deformation_type_t::quadratic)
        {
            object.integrate_quadratic(sim_params.dt);
        }
        else
        {
            object.integrate(sim_params.dt);
        }

        viewer.data().V = object.V();

        new_physics_update = true;
    };

    auto start = std::chrono::high_resolution_clock::now();

    viewer.callback_pre_draw =
        [&object_loaded, &new_physics_update, &simulate, &sim_params, &pick, &object, &start](
            igl::opengl::glfw::Viewer& viewer) -> bool {
        auto const now = std::chrono::high_resolution_clock::now();

        if (!object_loaded)
            return false;

        if (!sim_params.pause)
        {
            simulate();
        }

        if (pick.picked)
        {
            Eigen::Vector3d p1{
                pick.object->V()(pick.vi, 0),
                pick.object->V()(pick.vi, 1),
                pick.object->V()(pick.vi, 2)};

            double const x2 = static_cast<double>(viewer.current_mouse_x);
            double const y2 =
                viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);

            Eigen::Vector3d const p2 = igl::unproject(
                                           Eigen::Vector3f(x2, y2, 1.0f),
                                           viewer.core().view,
                                           viewer.core().proj,
                                           viewer.core().viewport)
                                           .cast<double>();

            Eigen::Vector3d d = (p2 - p1).normalized() * sim_params.PickAmplifier;
            pick.object->V().block(pick.vi, 0, 1, 3) += d.transpose();
        }

        // draw fixed points
        viewer.data().clear_points();
        Eigen::MatrixX3d P(object.count_fixed(), 3);
        Eigen::MatrixX3d C(object.count_fixed(), 3);
        for (unsigned int i = 0, j = 0; i < object.V().rows(); ++i)
        {
            if (!object.is_fixed(i))
                continue;
            P.block(j, 0, 1, 3) = object.V().block(i, 0, 1, 3);
            C.block(j, 0, 1, 3) = Eigen::RowVector3d(255., 0., 0.);
            ++j;
        }
        viewer.data().add_points(P, C);

        if (new_physics_update)
        {
            viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
            new_physics_update = false;
        }

        viewer.data().show_lines = sim_params.show_triangles;

        auto const end = std::chrono::high_resolution_clock::now();
        auto const period_in_nanoseconds =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        start = end;

        sim_params.fps          = 1'000'000'000 / period_in_nanoseconds;
        sim_params.vertices     = object.V().rows();
        sim_params.triangles    = object.F().rows();
        sim_params.total_memory = object.get_memory_in_bytes();

        return false;
    };

    viewer.launch();
    return EXIT_SUCCESS;
}
