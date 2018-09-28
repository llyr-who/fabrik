
#ifndef GLNIXAP_H
#define GLNIXAP_H
#include<iostream>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include<vector>
#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/extensions/Xrender.h>
#include <X11/Xutil.h>
#include<X11/X.h>
#include<X11/Xlib.h>
#include<GL/gl.h>
#include<GL/glx.h>
#include<GL/glu.h>
#include"GeometryGen.h"

#include "Timer.h"

#define USE_CHOOSE_FBCONFIG



class GLnixAPP
{
public:
	GLnixAPP();
	//virtual ~GLnixAPP();
	
	Window   MainWnd()const;
	float     AspectRatio()const;
	
	int Run();
 
	// Framework methods.  Derived client class overrides these methods to 
	// implement specific application requirements.

	virtual bool Init();
	//virtual void OnResize(); 
	virtual void UpdateScene(float dt);
	virtual void RedrawTheWindow();
	virtual int  UpdateTheMessageQueue();

	// Convenience overrides for handling mouse input.
	virtual void OnMouseDown(XButtonEvent btn,int x, int y){ }
	virtual void OnMouseUp(XButtonEvent btn,int x, int y)  { }
	virtual void OnMouseMove(int x, int y){ }

protected:
	int InitMainWindow();
	int InitGLnix();

	//void CalculateFrameStats();

protected:

	bool      appPaused;
	bool      minimized;
	bool      maximized;
	bool      resizing;

	GameTimer timer;

	int Xscreen;
	Atom del_atom;
	Colormap cmap;
	Display *Xdisplay;
	XVisualInfo *visual;
	XRenderPictFormat *pict_format;
	GLXFBConfig *fbconfigs, fbconfig;
	int numfbconfigs;
	GLXContext render_context;
	Window Xroot, window_handle;
	GLXWindow glX_window_handle;
	

	// Derived class should set these in derived constructor to customize starting values.
	int width, height;

};

#endif 
