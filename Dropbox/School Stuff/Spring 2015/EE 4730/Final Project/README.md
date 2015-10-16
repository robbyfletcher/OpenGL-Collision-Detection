Compilation: 

Included is an XCode Project. If you have access to a mac with XCode, you should just be able to open it and all the settings should be included. 
ls

If not, included is a main.cpp file in the finalProject folder, which is the entire code, with the exception of the Meshlib_Core files. To compile correctly, you should link both the Mesh.cpp file and the Vertex.cpp file with the compiler. I would recommend gcc. 

Running:

You must pass a mesh file as an argument to the program in order for it to run properly. If not, it will not run. 

Once the program is running you should see two meshes. The controls for moving the camera have not changed, however, pressing 'm' will enter 'move' mode, which allows you to move the second mesh to simulate collisions. 