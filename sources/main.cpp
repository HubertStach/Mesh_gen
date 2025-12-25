#include "raylib.h"

#include "geometry.h"

#include <vector>

#define SCREEN_WIDTH (1280)
#define SCREEN_HEIGHT (720)

#define WINDOW_TITLE "MESH_GEN"

/*
First time build:
mkdir build
cd build
cmake ..
cmake --build .
./Debug/mesh_gen.exe

After that:
cmake --build .
./Debug/mesh_gen.exe
*/

int main()
{
    float spacing = 100;
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, WINDOW_TITLE);
    SetTargetFPS(60);

    geo::Mesh mesh;

    bool mesh_created = false;

    while (!WindowShouldClose())
    {
        BeginDrawing();
         ClearBackground(BLACK);

        if(IsMouseButtonPressed(MOUSE_BUTTON_LEFT)){
            if(mesh_created){
                mesh.edges.clear();
                mesh.triangles.clear();
                mesh.nodes.clear();
                mesh.initial_bc_nodes.clear();
                mesh_created = false;
            }

            float x = GetMousePosition().x;
            float y = GetMousePosition().y;
            mesh.add_point(x, y);
        }

        if(IsKeyPressed(KEY_ENTER)){
            mesh.create_mesh(spacing);
            mesh_created = true;
        }

        //rysowanie

        mesh.draw_nodes(2);
        mesh.draw_tr();


        EndDrawing();
    }

    CloseWindow();

    return 0;
}
