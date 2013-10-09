Project_D
=========

Project D.A.V.E. (Drawing ~~and Animating~~ Awesome Virtual Environments)

Building
=========

The lib directory contains pre-build dependencies for some compilers.
If there aren't any for yours, you'll need to build the deps.

Build each dep in a subfolder 'build' to avoid polluting the source tree.

- GLFW:
GLFW is built via cmake. On Windows, you'll want to specify a generator like
  -G "Visual Studio 10"

Building with a makefile
-------------------------

There are some compiler / system specific makefiles in the git root (they point to the right lib directory). Rename one to 'Makefile' or make a new one.

Building with Visual Studio
---------------------------

There is a VS2010 solution which _should_ work. If you want to import this into a different version of VS, MAKE A COPY and put it into an appropriate directory (i.e. don't go overwriting my VS2010 solution!). The gitignore will also need updating if you want to commit a new VS solution.

Otherwise: Make a solution (empty C++ project) in the git root, then rename the solution dir to (eg) 'VS2010'.
- Add $(SolutionDir)../include to the include path
- Add $(SolutionDir)../lib/[APPROPRIATE_DIRECTORY] to the lib path
- Add the dependency lib names to the linker (leave the 'inherit' option on)
- Set the working directory to be $(SolutionDir)..
- Add the contents of the src directory to the project
- Use the same CRT as the deps were built with (default Multi-threaded Debug DLL)


