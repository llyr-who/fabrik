#include "GLnixAPP.h"

#include <sstream>

static Bool WaitForMapNotify(Display *d, XEvent *e, char *arg) {
    return d && e && arg && (e->type == MapNotify) &&
           (e->xmap.window == *(Window *)arg);
}

static void describe_fbconfig(GLXFBConfig fbconfig, Display *disp) {
    int doublebuffer;
    int red_bits, green_bits, blue_bits, alpha_bits, depth_bits;

    glXGetFBConfigAttrib(disp, fbconfig, GLX_DOUBLEBUFFER, &doublebuffer);
    glXGetFBConfigAttrib(disp, fbconfig, GLX_RED_SIZE, &red_bits);
    glXGetFBConfigAttrib(disp, fbconfig, GLX_GREEN_SIZE, &green_bits);
    glXGetFBConfigAttrib(disp, fbconfig, GLX_BLUE_SIZE, &blue_bits);
    glXGetFBConfigAttrib(disp, fbconfig, GLX_ALPHA_SIZE, &alpha_bits);
    glXGetFBConfigAttrib(disp, fbconfig, GLX_DEPTH_SIZE, &depth_bits);

    // we need to write some code to describe this to the user.
}

static int isExtensionSupported(const char *extList, const char *extension) {
    const char *start;
    const char *where, *terminator;

    /* Extension names should not have spaces. */
    where = strchr(extension, ' ');
    if (where || *extension == '\0') return 0;

    /* It takes a bit of care to be fool-proof about parsing the
       OpenGL extensions string. Don't be fooled by sub-strings,
       etc. */
    for (start = extList;;) {
        where = strstr(start, extension);

        if (!where) break;

        terminator = where + strlen(extension);

        if (where == start || *(where - 1) == ' ')
            if (*terminator == ' ' || *terminator == '\0') return 1;

        start = terminator;
    }
    return 0;
}

GLnixAPP::GLnixAPP()
    : width(1920),
      height(1080),
      appPaused(false),
      minimized(false),
      maximized(false),
      resizing(false){};

Window GLnixAPP::MainWnd() const { return window_handle; }

float GLnixAPP::AspectRatio() const {
    return static_cast<float>(width) / height;
}

int GLnixAPP::Run() {
    timer.Reset();

    while (UpdateTheMessageQueue()) {
        timer.Tick();

        if (!appPaused) {
            // CalculateFrameStats();
            UpdateScene(timer.DeltaTime());
            RedrawTheWindow();
        } else {
            // Sleep(100);
        }
    }
    return 1;
}

bool GLnixAPP::Init() {
    if (!InitMainWindow()) return false;

    if (!InitGLnix()) return false;

    return true;
}

// void GLnixAPP::OnResize()
//{
//	// use this function to deal with veiwportsa etc
//}

void GLnixAPP::UpdateScene(float dt){};
void GLnixAPP::RedrawTheWindow(){};
int GLnixAPP::UpdateTheMessageQueue() {
    XEvent event;
    XConfigureEvent *xc;

    while (XPending(Xdisplay)) {
        XNextEvent(Xdisplay, &event);
        switch (event.type) {
            case ClientMessage:
                if (event.xclient.data.l[0] == del_atom) {
                    return 0;
                }
                break;

            case ConfigureNotify:
                xc = &(event.xconfigure);
                width = xc->width;
                height = xc->height;
                break;

            case KeyPress:
                break;
            case MotionNotify:
                OnMouseMove(event.xbutton.x, event.xbutton.y);
                break;
            case ButtonPress:
                OnMouseDown(event.xbutton, event.xbutton.x, event.xbutton.y);
                break;
            case ButtonRelease:
                OnMouseUp(event.xbutton, event.xbutton.x, event.xbutton.y);
                break;
        }
    }
    return 1;
}

int GLnixAPP::InitMainWindow() {
    static int VisData[] = {GLX_RENDER_TYPE,
                            GLX_RGBA_BIT,
                            GLX_DRAWABLE_TYPE,
                            GLX_WINDOW_BIT,
                            GLX_DOUBLEBUFFER,
                            True,
                            GLX_RED_SIZE,
                            8,
                            GLX_GREEN_SIZE,
                            8,
                            GLX_BLUE_SIZE,
                            8,
                            GLX_ALPHA_SIZE,
                            8,
                            GLX_DEPTH_SIZE,
                            16,
                            None};
    XEvent event;
    int x, y, attr_mask;
    XSizeHints hints;
    XWMHints *startup_state;
    XTextProperty textprop;
    XSetWindowAttributes attr = {
        0,
    };
    static char *title = "GLnix Engine";

    Xdisplay = XOpenDisplay(NULL);
    if (!Xdisplay) {
        std::cout << "Couldn't connect to X server\n";
    }
    Xscreen = DefaultScreen(Xdisplay);
    Xroot = RootWindow(Xdisplay, Xscreen);

    fbconfigs = glXChooseFBConfig(Xdisplay, Xscreen, VisData, &numfbconfigs);

    fbconfig = 0;
    for (int i = 0; i < numfbconfigs; i++) {
        visual =
            (XVisualInfo *)glXGetVisualFromFBConfig(Xdisplay, fbconfigs[i]);
        if (!visual) continue;

        pict_format = XRenderFindVisualFormat(Xdisplay, visual->visual);
        if (!pict_format) continue;

        fbconfig = fbconfigs[i];
        if (pict_format->direct.alphaMask > 0) {
            break;
        }
    }

    if (!fbconfig) {
        std::cout << "No matching FB config found";
    }

    describe_fbconfig(fbconfig, Xdisplay);

    cmap = XCreateColormap(Xdisplay, Xroot, visual->visual, AllocNone);

    attr.colormap = cmap;
    attr.background_pixmap = None;
    attr.border_pixmap = None;
    attr.border_pixel = 0;
    attr.event_mask = StructureNotifyMask | EnterWindowMask | LeaveWindowMask |
                      ExposureMask | ButtonPressMask | ButtonReleaseMask |
                      OwnerGrabButtonMask | KeyPressMask | KeyReleaseMask |
                      ButtonMotionMask;

    attr_mask =
        //	CWBackPixmap|
        CWColormap | CWBorderPixel | CWEventMask;

    x = width / 2, y = height / 2;

    window_handle =
        XCreateWindow(Xdisplay, Xroot, x, y, width, height, 0, visual->depth,
                      InputOutput, visual->visual, attr_mask, &attr);

    if (!window_handle) {
        std::cout << "Couldn't create the window\n";
    }

    glX_window_handle = window_handle;

    textprop.value = (unsigned char *)title;
    textprop.encoding = XA_STRING;
    textprop.format = 8;
    textprop.nitems = strlen(title);

    hints.x = x;
    hints.y = y;
    hints.width = width;
    hints.height = height;
    hints.flags = USPosition | USSize;

    startup_state = XAllocWMHints();
    startup_state->initial_state = NormalState;
    startup_state->flags = StateHint;

    XSetWMProperties(Xdisplay, window_handle, &textprop, &textprop, NULL, 0,
                     &hints, startup_state, NULL);

    XFree(startup_state);

    XMapWindow(Xdisplay, window_handle);
    XIfEvent(Xdisplay, &event, WaitForMapNotify, (char *)&window_handle);

    if ((del_atom = XInternAtom(Xdisplay, "WM_DELETE_WINDOW", 0)) != None) {
        XSetWMProtocols(Xdisplay, window_handle, &del_atom, 1);
    }
    return 1;
}

int GLnixAPP::InitGLnix() {
    // nice trick here
    int dummy;
    if (!glXQueryExtension(Xdisplay, &dummy, &dummy)) {
        std::cout << "OpenGL not supported by X server\n";
    }

#if USE_GLX_CREATE_CONTEXT_ATTRIB
#define GLX_CONTEXT_MAJOR_VERSION_ARB 0x2091
#define GLX_CONTEXT_MINOR_VERSION_ARB 0x2092
    render_context = NULL;
    if (isExtensionSupported(
            glXQueryExtensionsString(Xdisplay, DefaultScreen(Xdisplay)),
            "GLX_ARB_create_context")) {
        typedef GLXContext (*glXCreateContextAttribsARBProc)(
            Display *, GLXFBConfig, GLXContext, Bool, const int *);
        glXCreateContextAttribsARBProc glXCreateContextAttribsARB =
            (glXCreateContextAttribsARBProc)glXGetProcAddressARB(
                (const GLubyte *)"glXCreateContextAttribsARB");
        if (glXCreateContextAttribsARB) {
            int context_attribs[] = {GLX_CONTEXT_MAJOR_VERSION_ARB, 3,
                                     GLX_CONTEXT_MINOR_VERSION_ARB, 0,
                                     // GLX_CONTEXT_FLAGS_ARB        ,
                                     // GLX_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB,
                                     None};

            int (*oldHandler)(Display *, XErrorEvent *) =
                XSetErrorHandler(&ctxErrorHandler);

            render_context = glXCreateContextAttribsARB(Xdisplay, fbconfig, 0,
                                                        True, context_attribs);

            XSync(Xdisplay, False);
            XSetErrorHandler(oldHandler);

            fputs("glXCreateContextAttribsARB failed", stderr);
        } else {
            fputs("glXCreateContextAttribsARB could not be retrieved", stderr);
        }
    } else {
        fputs("glXCreateContextAttribsARB not supported", stderr);
    }

    if (!render_context) {
#else
    {
#endif
        render_context =
            glXCreateNewContext(Xdisplay, fbconfig, GLX_RGBA_TYPE, 0, True);
        if (!render_context) {
            std::cout << "Failed to create a GL context\n";
        }
    }

    if (!glXMakeContextCurrent(Xdisplay, glX_window_handle, glX_window_handle,
                               render_context)) {
        std::cout << "glXMakeCurrent failed for window\n";
    }

    return 1;
}

