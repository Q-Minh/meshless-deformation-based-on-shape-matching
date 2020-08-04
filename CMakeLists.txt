cmake_minimum_required(VERSION 3.14)
project(MeshlessDeformationsBasedOnShapeMatching)

add_executable(App)
set_target_properties(App PROPERTIES CXX_STANDARD 14)

find_package(Eigen3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
target_link_libraries(App PRIVATE glfw glad::glad Eigen3::Eigen)

target_sources(App 
PRIVATE
# headers
	src/deformable_mesh.h

	external/imgui/imgui_fonts_droid_sans.h
	external/imgui/imgui_impl_glfw.h
	external/imgui/imgui_impl_opengl3.h

	external/imgui/imconfig.h
	external/imgui/imgui.cpp
	external/imgui/imgui.h
	external/imgui/imgui_demo.cpp
	external/imgui/imgui_draw.cpp
	external/imgui/imgui_internal.h
	external/imgui/imgui_widgets.cpp
	external/imgui/imstb_rectpack.h
	external/imgui/imstb_textedit.h
	external/imgui/imstb_truetype.h

	external/imgui/misc/cpp/imgui_stdlib.cpp
	external/imgui/misc/cpp/imgui_stdlib.h

# sources
    src/main.cpp
	src/deformable_mesh.cpp

	external/imgui/imgui_impl_glfw.cpp
	external/imgui/imgui_impl_opengl3.cpp
)

target_include_directories(App
PRIVATE
	src/
	external/
	external/imgui/
	libigl/include/
)

target_compile_options(App PUBLIC /bigobj)