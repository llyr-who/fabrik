#pragma once

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glxext.h>

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

namespace ext {
void load() {
    const unsigned char* glssp =
        reinterpret_cast<const unsigned char*>("glShaderSource");
    glShaderSource = (PFNGLSHADERSOURCEPROC)glXGetProcAddress(glssp);

    const unsigned char* glcmplsp =
        reinterpret_cast<const unsigned char*>("glCompileShader");
    glCompileShader = (PFNGLCOMPILESHADERPROC)glXGetProcAddress(glcmplsp);

    const unsigned char* glcsp =
        reinterpret_cast<const unsigned char*>("glCreateShader");
    glCreateShader = (PFNGLCREATESHADERPROC)glXGetProcAddress(glcsp);

    const unsigned char* glasp =
        reinterpret_cast<const unsigned char*>("glAttachShader");
    glAttachShader = (PFNGLATTACHSHADERPROC)glXGetProcAddress(glasp);

    const unsigned char* glcspp =
        reinterpret_cast<const unsigned char*>("glCreateProgram");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC)glXGetProcAddress(glcsp);

    const unsigned char* glevaap =
        reinterpret_cast<const unsigned char*>("glEnableVertexAttribArray");
    glEnableVertexAttribArray =
        (PFNGLENABLEVERTEXATTRIBARRAYPROC)glXGetProcAddress(glevaap);

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
}  // namespace ext
