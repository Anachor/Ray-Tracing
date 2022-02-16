# Ray Tracing Project

This project is a ray-tracing engine designed to render 3D scenes with realistic lighting, shading, and reflections. It supports various geometric objects, recursive ray tracing, and customizable material properties.

## Features

- **Objects**:
  - Spheres
  - Triangles
  - General Quadric Surfaces
  - Floor with alternating tile colors
- **Lighting**:
  - Ambient, diffuse, and specular lighting
  - Multiple light sources
- **Materials**:
  - Configurable reflection coefficients (ambient, diffuse, specular, recursive)
  - Adjustable shininess for specular highlights
- **Camera**:
  - Configurable camera positions and orientations
- **Rendering**:
  - Recursive ray tracing for reflections
  - Adjustable recursion depth
  - Pixel-based rendering with configurable resolution

## File Descriptions

### `scene.txt`
Defines the scene to be rendered, including:
- Number of objects and light sources
- Object properties (position, size, color, material coefficients)
- Light source properties (position, color)
- Camera configuration and rendering parameters

### `objects.h`
Contains the implementation of:
- Geometric objects (spheres, triangles, quadric surfaces, floor)
- Lighting and shading calculations
- Ray-object intersection logic
- Recursive ray tracing algorithm

### `geometry.h`
Provides geometric utilities for vector and point operations, such as:
- Dot and cross products
- Distance calculations
- Plane and line equations

## How to Run

1. Define the scene in `scene.txt`.
2. Compile the project using a `g++ -O2 main.cpp -lglut -lGL -lGLU`
3. Run the executable to render the scene.



## Acknowledgments

This project is part of a ray-tracing assignment and uses OpenGL for rendering. Additionally, the project makes use of the `bitmap_image` library for image generation and output.