#include"GLnixAPP.h"
#include"myextloader.h"
#include"GeometryGen.h"
#include"MathHelper.h"
#include"cloth.h"
using namespace std;
float const light0_dir[]={0,1,1,0};
float const light0_color[]={78./255., 8./255., 184./255.,1};

float const light1_dir[]={-1,1,1,0};
float const light1_color[]={25./255., 220./255., 70./255.,1};

float const light2_dir[]={0,-1,0,0};
float const light2_color[]={31./255., 75./255., 160./255.,1};


class GLnixDEMO : public GLnixAPP
{
public:
	GLnixDEMO();


	bool Init();
	void UpdateScene(float dt);
	void RedrawTheWindow(); 

	void OnMouseDown(XButtonEvent btn,int x, int y);
	void OnMouseUp(XButtonEvent btn,int x, int y);
	void OnMouseMove(int x, int y);


private:
	void BuildClothGeometryBuffers();
private:
	UINT sphereVertexOffset;
	UINT sphereIndexOffset;
	UINT sphereIndexCount;
	
	UINT VAOcloth;

	UINT BUFFERS[2];
	
	Cloth cloth;
	UINT gridIndexCount;
	int mousex;
	int mousey;
	int but;
	
	AV4X4FLOAT viewModelMatrix;
	AV4X4FLOAT projMatrix;

	AV4X4FLOAT gridWorldMatrix;

	float theta;
	float phi;
	float radius;

};


int main()
{   
	load_extension_function_pointers();

	GLnixDEMO theApp;
	if( !theApp.Init() )
		return 0;
	return theApp.Run();
}

GLnixDEMO::GLnixDEMO()
: GLnixAPP(), phi(1.5f*MathHelper::Pi),theta(1.5f*MathHelper::Pi), radius(10)
{
	mousex = 0;
	mousey = 0;
	AV4FLOAT r(1,1,1,1);
	AV4X4FLOAT I;
	I.diag (r);
 
	gridWorldMatrix = I;
	projMatrix = I;
	viewModelMatrix = I;

}



bool GLnixDEMO::Init()
{
	

	if(!GLnixAPP::Init())
		return false;

	//cloth.Init(200,200,0.02,0.01,1000,1200,0.1,0.4,0.4);
	//cloth.Init(300,300,0.06,0.007,1000,1200,0.1,0.4,0.4);
	//cloth.Init(200,200,0.03,0.01,2000,1800,0.1,0.3,0.6); this one is really good  A make red for "spawns cape"
	//cloth.Init(200,200,0.03,0.01,2000,1800,0.2,0.3,0.8); this one is also really good B1 .. this one is the same deal.
	// cloth.Init(200,200,0.025,0.01,2500,2000,1.0,2.0,0.9); C make this one white for a bedsheet look, (works well with gravity too) 
	//cloth.Init(200,200,0.025,0.01,3000,2500,1.0,2.0,0.9); variation of C

	cloth.Init(200,200,0.025,0.01,3000,2500,1.1,1.5,0.8); //Blue velvet, I like it
	BuildClothGeometryBuffers();
	projMatrix =formProjMatrix(0.25f*MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	return true;
}



void GLnixDEMO::UpdateScene(float dt)
{

	float x = radius*sinf(phi)*cosf(theta);
	float z = radius*sinf(phi)*sinf(theta);
	float y = radius*cosf(phi);


	AV4FLOAT position(x,y,z,1.0);
	
	AV4FLOAT target(0.0,0.0,0.0,0.0);
	AV4FLOAT up(0.0,1.0,0.0,0.0);

	viewModelMatrix = formViewModelMatrix(position,target,up);

    float T = timer.GameTime();
	
/*	cloth.Update(dt,2*cos(T)*cos(T),cos(T)*cos(T),0); when grav engaged???*/
	//cloth.Update(dt,2*cos(T)*cos(T),0,2*sin(T)*cos(T));

	
	//cloth.Update(dt,MathHelper::RandF(0.9,1.0),2*cos(3*T),MathHelper::RandF(0.1,0.2)); A 


	cloth.Update(dt,MathHelper::RandF(0.9,1.0),MathHelper::RandF(1.0,3.0)*cos(3*T),MathHelper::RandF(0.1,0.2)); // C/B

	
   	GLnix_glBindVertexArray(VAOcloth);
	vector<GeometryGenerator::Vertex> v(cloth.VertexCount());
	for(int i =0; i < cloth.VertexCount();i++)
	{
		v[i].Position = cloth[i];
		v[i].Normal = cloth.Normal(i);
	}
	GLnix_glBufferSubData(GL_ARRAY_BUFFER,
	                      0,
	                      cloth.VertexCount() * sizeof(GLfloat) * 11,
	                      &v[0]);

}



void GLnixDEMO::RedrawTheWindow()
{
	float const aspect = AspectRatio();

	float x = radius*sinf(phi)*cosf(theta);
	float z = radius*sinf(phi)*sinf(theta);
	float y = radius*cosf(phi);



	glDrawBuffer(GL_BACK);

	glViewport(0, 0, width, height);
	glClearColor(0.0, 0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glLoadMatrixf(projMatrix.m);

	


	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLightfv(GL_LIGHT0, GL_POSITION, light0_dir);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_color);

	glLightfv(GL_LIGHT1, GL_POSITION, light1_dir);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_color);

	glLightfv(GL_LIGHT2, GL_POSITION, light2_dir);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_color);


	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glCullFace(GL_BACK);

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(viewModelMatrix.m);
	glColor4f(0.1, 0.2, 0.8, 1);
	GLnix_glBindVertexArray(VAOcloth);
    glDrawElements(GL_TRIANGLES,  3*cloth.TriangleCount(), GL_UNSIGNED_INT, 0);
    glCullFace(GL_FRONT);
	GLnix_glBindVertexArray(VAOcloth);
    glDrawElements(GL_TRIANGLES,  3*cloth.TriangleCount(), GL_UNSIGNED_INT, 0);
 	glXSwapBuffers(Xdisplay, glX_window_handle);



}




void GLnixDEMO::BuildClothGeometryBuffers()
{


	GeometryGenerator::MeshData grid;
	GeometryGenerator geoGen;
	geoGen.CreateGrid(200.0f, 200.0f, 50, 50, grid);

	std::vector<UINT> indices(3*cloth.TriangleCount()); // 3 indices per face

	// Iterate over each quad.
	UINT m = cloth.RowCount();
	UINT n = cloth.ColumnCount();
	int k = 0;
	for(UINT i = 0; i < m-1; ++i)
	{
		for(UINT j = 0; j < n-1; ++j)
		{
			indices[k]   = i*n+j;
			indices[k+1] = i*n+j+1;
			indices[k+2] = (i+1)*n+j;

			indices[k+3] = (i+1)*n+j;
			indices[k+4] = i*n+j+1;
			indices[k+5] = (i+1)*n+j+1;

			k += 6; // next quad
		}
	}


	GLnix_glGenVertexArrays(1,&VAOcloth);
	GLnix_glGenBuffers(2,BUFFERS);
	GLnix_glBindVertexArray(VAOcloth); 
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);


	GLnix_glBindBuffer(GL_ARRAY_BUFFER,BUFFERS[0]);
	GLnix_glBufferData(GL_ARRAY_BUFFER,
	                              cloth.VertexCount() * sizeof(GLfloat) * 11,
	                              0, GL_STATIC_DRAW);


	GLnix_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BUFFERS[1]);
    GLnix_glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                                   indices.size() * sizeof(UINT), 
                                   &indices.front(), GL_STATIC_DRAW);
	
    glVertexPointer(3, GL_FLOAT,sizeof(GLfloat) * 11, 0);
	glNormalPointer(GL_FLOAT, sizeof(GLfloat) * 11,  (GLvoid*)(3*sizeof(GLfloat)) );
	


}






void GLnixDEMO::OnMouseDown(XButtonEvent btn,int x, int y)
{
	mousex = x;
	mousey = y;
	but = btn.button;
}
void GLnixDEMO::OnMouseUp(XButtonEvent btn,int x, int y)
{

}
void GLnixDEMO::OnMouseMove(int x, int y)
{
	if(but == 1)
	{
		// Make each pixel correspond to a quarter of a degree.
		float dx = ANTMATHConvertToRadians(0.25f*static_cast<float>(x - mousex));
		float dy = ANTMATHConvertToRadians(0.25f*static_cast<float>(y -mousey));

		// Update angles based on input to orbit camera around box.
		theta += dx;
		phi   += dy;

		// Restrict the angle mPhi.
		phi = MathHelper::Clamp(phi, 0.1f, MathHelper::Pi-0.1f);
	}
	else if (but == 3)
	{
		// Make each pixel correspond to 0.2 unit in the scene.
		float dx = 0.2f*static_cast<float>(x - mousex);
		float dy = 0.2f*static_cast<float>(y - mousey);

		// Update the camera radius based on input.
		radius += dx - dy;

		// Restrict the radius.
		radius = MathHelper::Clamp(radius, 1.0f, 20.0f);
	}

	

	mousex = x;
	mousey = y;
}

