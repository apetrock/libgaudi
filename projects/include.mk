GLVUFLAGS = -msse3 -lglvuglvu -lglvucamera -lglvufileutils -lglvumath -L${GLVU_LIB}

#bullet opengl -resolve later
#LDFLAGS_COMMON =  -lGLEW -lGL -lGLU  -lglut -lX11 -L/usr/local/lib/ -L/usr/local/lib/x86_64 -L/opt/local/lib  -lode -lbulletopenglsupport -L../../src -lBulletDynamics -lBulletCollision -lLinearMath -L/home/john/Code/libgaudi/src/ode-0.13/ode/build/ -lstdc++ -lpng -lz -ljpeg -fopenmp -msse2

#-lnanogui

PKG_LIBS = $(shell pkg-config --static --libs x11 xrandr xinerama xcursor xi xxf86vm glew glfw3)

ARPACK_INC = -I../../../arpackpp/include -I../../../arpackpp/examples/matrices/nonsym/ -I../../../arpackpp/examples/matrices/sym/ -I../../../arpackpp/examples/areig/  -I../../../arpackpp/examples/areig/sym/ -I../../../arpackpp/external/SuperLU/SRC

ARPACK_LIBS = -L/usr/lib/ -llapack -L../../../arpackpp/external/ -lopenblas -lsuperlu -larpack 

LDFLAGS_COMMON =  -lGL -lGLU $(PKG_LIBS) $(ARPACK_LIBS)  -lpthread -ldl -lnanogui -L../../src/nanogui/src   -L/user/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -L/usr/local/lib/x86_64 -L/opt/local/lib -L../../src  -lstdc++ -lpng -lz -ljpeg -fopenmp -msse2

#-I../../src/Eigen/ -I../../src/volume/ -I../../src/gaudi/ -I../../src/ode-0.13/include -I../../src/ode-0.13/include/drawstuff

CFLAGS_COMMON = -c -Wall -DDO_PNG_OUT=0 -I./  -I../../ -I../../src -I../../src/nanogui/include -I../../src/nanoguihelpers -I../../src/eigen -I../../src/gaudi/ -I../../src/quartic/ -I../../src/zither/ -I/opt/local/include -I/usr/local/include/Eigen  -O3 -DNO_FFT -fopenmp -msse2

#NANO_SRC = ../../src/nanogui/nanovg/nanovg.c \
	../../src/nanogui/resources.cpp \
	../../src/nanogui/src/nanogui.cpp \
	../../src/nanogui/src/glutil.cpp \
	../../src/nanogui/src/screen.cpp \
	../../src/nanogui/src/window.cpp \
	../../src/nanogui/src/widget.cpp \
	../../src/nanogui/src/layout.cpp \
	../../src/nanogui/src/label.cpp \
	../../src/nanogui/src/theme.cpp \
	../../src/nanogui/src/checkbox.cpp \
	../../src/nanogui/src/button.cpp \
	../../src/nanogui/src/popup.cpp \
	../../src/nanogui/src/popupbutton.cpp \
	../../src/nanogui/src/combobox.cpp \
	../../src/nanogui/src/progressbar.cpp \
	../../src/nanogui/src/messagedialog.cpp \
	../../src/nanogui/src/textbox.cpp \
	../../src/nanogui/src/slider.cpp \
	../../src/nanogui/src/imagepanel.cpp \
	../../src/nanogui/src/imageview.cpp \
	../../src/nanogui/src/vscrollpanel.cpp \
	../../src/nanogui/src/colorwheel.cpp \
