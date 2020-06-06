#pragma once

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glxext.h>

extern PFNGLBINDBUFFERPROC glBindBuffer;
extern PFNGLCREATESHADERPROC glCreateShader;
extern PFNGLSHADERSOURCEPROC glShaderSource;
extern PFNGLGENBUFFERSPROC glGenBuffers;
extern PFNGLDELETEBUFFERSPROC glDeleteBuffer;
extern PFNGLBUFFERDATAPROC glBufferData;
extern PFNGLBUFFERSUBDATAPROC glBufferSubData;
extern PFNGLMAPBUFFERPROC glMapBuffer;
extern PFNGLUNMAPBUFFERPROC glUnmapBuffer;
extern PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
extern PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;

void load_extension_function_pointers() {
    const unsigned char* glssp =
        reinterpret_cast<const unsigned char*>("glShaderSource");
    glShaderSource = (PFNGLSHADERSOURCEPROC)glXGetProcAddress(glssp);

    const unsigned char* glcsp =
        reinterpret_cast<const unsigned char*>("glCreateShader");
    glCreateShader = (PFNGLCREATESHADERPROC)glXGetProcAddress(glcsp);

    const unsigned char* tt =
        reinterpret_cast<const unsigned char*>("glBindBuffer");
    glBindBuffer = (PFNGLBINDBUFFERPROC)glXGetProcAddress(tt);

    const unsigned char* tr =
        reinterpret_cast<const unsigned char*>("glGenBuffers");
    glGenBuffers = (PFNGLGENBUFFERSPROC)glXGetProcAddress(tr);

    const unsigned char* tf =
        reinterpret_cast<const unsigned char*>("glDeleteBuffer");
    glDeleteBuffer = (PFNGLDELETEBUFFERSPROC)glXGetProcAddress(tf);

    const unsigned char* ff =
        reinterpret_cast<const unsigned char*>("glBufferData");
    glBufferData = (PFNGLBUFFERDATAPROC)glXGetProcAddress(ff);

    const unsigned char* mm =
        reinterpret_cast<const unsigned char*>("glMapBuffer");
    glMapBuffer = (PFNGLMAPBUFFERPROC)glXGetProcAddress(mm);

    const unsigned char* ml =
        reinterpret_cast<const unsigned char*>("glUnmapBuffer");
    glUnmapBuffer = (PFNGLUNMAPBUFFERPROC)glXGetProcAddress(ml);

    const unsigned char* mva =
        reinterpret_cast<const unsigned char*>("glBindVertexArray");
    glBindVertexArray = (PFNGLBINDVERTEXARRAYPROC)glXGetProcAddress(mva);

    const unsigned char* mvb =
        reinterpret_cast<const unsigned char*>("glGenVertexArrays");
    glGenVertexArrays = (PFNGLGENVERTEXARRAYSPROC)glXGetProcAddress(mvb);

    const unsigned char* mvbs =
        reinterpret_cast<const unsigned char*>("glBufferSubData");
    glBufferSubData = (PFNGLBUFFERSUBDATAPROC)glXGetProcAddress(mvbs);
}
