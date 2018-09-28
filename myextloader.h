
#ifndef MYEXTENSIONLOADER_H
#define MYEXTENSIONLOADER_H

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glxext.h>

PFNGLBINDBUFFERPROC GLnix_glBindBuffer;
PFNGLGENBUFFERSPROC GLnix_glGenBuffers;
PFNGLDELETEBUFFERSPROC GLnix_glDeleteBuffer;
PFNGLBUFFERDATAPROC GLnix_glBufferData;
PFNGLBUFFERSUBDATAPROC GLnix_glBufferSubData;
PFNGLMAPBUFFERPROC GLnix_glMapBuffer;
PFNGLUNMAPBUFFERPROC GLnix_glUnmapBuffer;
PFNGLBINDVERTEXARRAYPROC GLnix_glBindVertexArray;
PFNGLGENVERTEXARRAYSPROC GLnix_glGenVertexArrays;

void load_extension_function_pointers()
{
 const unsigned char* tt = reinterpret_cast<const unsigned char *>( "glBindBuffer" );
 GLnix_glBindBuffer = (PFNGLBINDBUFFERPROC) glXGetProcAddress(tt);
	
const unsigned char* tr = reinterpret_cast<const unsigned char *>( "glGenBuffers" );
 GLnix_glGenBuffers = (PFNGLGENBUFFERSPROC) glXGetProcAddress(tr);
	
 const unsigned char* tf = reinterpret_cast<const unsigned char *>( "glDeleteBuffer" );
 GLnix_glDeleteBuffer = (PFNGLDELETEBUFFERSPROC) glXGetProcAddress(tf);

 const unsigned char* ff = reinterpret_cast<const unsigned char *>( "glBufferData" );
 GLnix_glBufferData = (PFNGLBUFFERDATAPROC) glXGetProcAddress(ff);

 const unsigned char* mm = reinterpret_cast<const unsigned char *>( "glMapBuffer" );
 GLnix_glMapBuffer = (PFNGLMAPBUFFERPROC) glXGetProcAddress(mm);


 const unsigned char* ml = reinterpret_cast<const unsigned char *>( "glUnmapBuffer" );
 GLnix_glUnmapBuffer = (PFNGLUNMAPBUFFERPROC) glXGetProcAddress(ml);
	
 const unsigned char* mva = reinterpret_cast<const unsigned char *>( "glBindVertexArray" );
 GLnix_glBindVertexArray = (PFNGLBINDVERTEXARRAYPROC) glXGetProcAddress(mva);

 const unsigned char* mvb = reinterpret_cast<const unsigned char *>( "glGenVertexArrays" );
 GLnix_glGenVertexArrays = (PFNGLGENVERTEXARRAYSPROC) glXGetProcAddress(mvb);


const unsigned char* mvbs = reinterpret_cast<const unsigned char *>( "glBufferSubData" );
 GLnix_glBufferSubData = (PFNGLBUFFERSUBDATAPROC) glXGetProcAddress(mvbs);
}






#endif
