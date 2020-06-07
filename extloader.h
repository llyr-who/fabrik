#pragma once

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glxext.h>

#include <stdexcept>

extern PFNGLBINDBUFFERPROC glBindBuffer;
extern PFNGLCREATESHADERPROC glCreateShader;
extern PFNGLATTACHSHADERPROC glAttachShader;
extern PFNGLLINKPROGRAMPROC glLinkProgram;
extern PFNGLUSEPROGRAMPROC glUseProgram;
extern PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer;
extern PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv;
extern PFNGLGETATTRIBLOCATIONPROC glGetAttribLocation;
extern PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
extern PFNGLCREATEPROGRAMPROC glCreateProgram;
extern PFNGLCOMPILESHADERPROC glCompileShader;
extern PFNGLSHADERSOURCEPROC glShaderSource;
extern PFNGLGENBUFFERSPROC glGenBuffers;
extern PFNGLDELETEBUFFERSPROC glDeleteBuffers;
extern PFNGLBUFFERDATAPROC glBufferData;
extern PFNGLBUFFERSUBDATAPROC glBufferSubData;
extern PFNGLMAPBUFFERPROC glMapBuffer;
extern PFNGLUNMAPBUFFERPROC glUnmapBuffer;
extern PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
extern PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;

void load() {
    const unsigned char* glssp =
        reinterpret_cast<const unsigned char*>("glShaderSource");
    glShaderSource = (PFNGLSHADERSOURCEPROC)glXGetProcAddress(glssp);
    if (!glShaderSource) {
        throw std::runtime_error("glShaderSource is nullptr");
    }

    const unsigned char* glcmplsp =
        reinterpret_cast<const unsigned char*>("glCompileShader");
    glCompileShader = (PFNGLCOMPILESHADERPROC)glXGetProcAddress(glcmplsp);
    if (!glCompileShader) {
        throw std::runtime_error("glCompile is nullptr");
    }

    const unsigned char* glcsp =
        reinterpret_cast<const unsigned char*>("glCreateShader");
    glCreateShader = (PFNGLCREATESHADERPROC)glXGetProcAddress(glcsp);

    const unsigned char* glupp =
        reinterpret_cast<const unsigned char*>("gtlUseProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC)glXGetProcAddress(glupp);

    const unsigned char* glasp =
        reinterpret_cast<const unsigned char*>("glAttachShader");
    glAttachShader = (PFNGLATTACHSHADERPROC)glXGetProcAddress(glasp);
    if (!glAttachShader) {
        throw std::runtime_error("glAttachShader is nullptr");
    }

    const unsigned char* gllpp =
        reinterpret_cast<const unsigned char*>("glLinkProgram");
    glLinkProgram = (PFNGLLINKPROGRAMPROC)glXGetProcAddress(gllpp);
    if (!glLinkProgram) {
        throw std::runtime_error("glLinkProgram is nullptr");
    }

    const unsigned char* glcspp =
        reinterpret_cast<const unsigned char*>("glCreateProgram");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC)glXGetProcAddress(glcspp);

    const unsigned char* glevaap =
        reinterpret_cast<const unsigned char*>("glEnableVertexAttribArray");
    glEnableVertexAttribArray =
        (PFNGLENABLEVERTEXATTRIBARRAYPROC)glXGetProcAddress(glevaap);

    const unsigned char* glgulp =
        reinterpret_cast<const unsigned char*>("glGetUniformLocation");
    glGetUniformLocation =
        (PFNGLGETUNIFORMLOCATIONPROC)glXGetProcAddress(glgulp);

    const unsigned char* glgalp =
        reinterpret_cast<const unsigned char*>("glGetAttribLocation");
    glGetAttribLocation = (PFNGLGETATTRIBLOCATIONPROC)glXGetProcAddress(glgalp);

    const unsigned char* glvapp =
        reinterpret_cast<const unsigned char*>("glVertexAttribPointer");
    glVertexAttribPointer =
        (PFNGLVERTEXATTRIBPOINTERPROC)glXGetProcAddress(glvapp);

    const unsigned char* glum4fvp =
        reinterpret_cast<const unsigned char*>("glUniformMatrix4fv");
    glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC)glXGetProcAddress(glum4fvp);

    const unsigned char* tt =
        reinterpret_cast<const unsigned char*>("glBindBuffer");
    glBindBuffer = (PFNGLBINDBUFFERPROC)glXGetProcAddress(tt);

    const unsigned char* tr =
        reinterpret_cast<const unsigned char*>("glGenBuffers");
    glGenBuffers = (PFNGLGENBUFFERSPROC)glXGetProcAddress(tr);

    const unsigned char* tf =
        reinterpret_cast<const unsigned char*>("glDeleteBuffers");
    glDeleteBuffers = (PFNGLDELETEBUFFERSPROC)glXGetProcAddress(tf);

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
