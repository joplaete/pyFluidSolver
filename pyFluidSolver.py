'''
    Created on Jul 8, 2010
    
    @author: jo plaete
    '''

import sys, math
try: from PySide import QtGui, QtCore
except: from PyQt4 import QtGui, QtCore
try: from scipy import weave
except: print "scipy.weave import failed - FluidSolverC won't work, pure python solver will be used"



class DrawFluidQt(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        
        self.dimension = 400
        self.n = 50
        self.rectSize = self.dimension / self.n
        self.dt = .3 # timestep: the higher the faster the fluid will evolve        
        
        # py or c extended solver
        try:
            from scipy import weave
            self.enableCExtentions(True)
        except:
            self.enableCExtentions(False)
        
        self.fs.setup(self.n, self.dt, viscosity=0.0, diffusion=0.0, vorticityConfinement=1, linearSolverIterations=5)
        # draw options
        self.drawVelocityField = False
        self.add_type = "density"
        # mouse position / cell index
        self.x = self.rectSize * self.n
        self.y = self.rectSize * self.n
        self.i = 0
        self.j = 0
        # qt
        self.setGeometry(300, 300, self.dimension, self.dimension)
        self.setWindowTitle('PyFluidSolver 1.0')
        # timer and solver update in other thread
        self.solverThread = SolverThread(self.fs)
        self.connect(self.solverThread, QtCore.SIGNAL('solverUpdate()'), self.repaint)
        self.solverThread.start()
        # image
        self.saveImage = False # on-off flag (shortcut i), note that images will always be written out at (non-scaled) true size where 1n=1pixel
        self.image = QtGui.QImage(self.dimension,self.dimension,QtGui.QImage.Format_RGB32)
        self.imagePath = "/tmp/"
        self.imageName = "pyFluid"
        self.imageType = "JPG" # JPG PNG GIF TIFF BMP
        self.frame = 0
        self.updateStatusBar()
    
    def enableCExtentions(self, enable):
        if enable:
            self.fs = FluidSolverC()
        else:
            self.fs = FluidSolver()
    
    def I(self,i,j):
        return i + (self.n + 2) * j
    
    def timerEvent(self, event):
        # update solver
        self.fs.velocitySolver()
        self.fs.densitySolver()
        self.repaint()
    
    def paintEvent_rects(self, event):
        """ Drawing rect for each cells, getting slower in high cell counts
            Can also draw velocity
            """
        paint = QtGui.QPainter()
        paint.begin(self)
        # draw grid
        rectSize = self.rectSize
        color = QtGui.QColor(255, 255, 255)
        paint.setPen(color)
        paint.setBrush(color)
        paint.drawRect(self.rectSize/2, self.rectSize/2, (self.n-1)*self.rectSize, (self.n-1)*self.rectSize)
        if self.saveImage: self.image.fill(color.rgb())
        for i in xrange(self.n):
            for j in xrange(self.n):
                dx = (i * (self.rectSize))
                dy = (j * (self.rectSize))
                # density
                if i == 0 or i == self.n-1 or j == 0 or j == self.n-1: continue
                elif self.fs.d[self.I(i, j)] > 0.0085: # making this a small value rather than 0 to avoid drawing those cells
                    # get denisty value from solver density array
                    c = int( (1.0 - self.fs.d[self.I(i, j)]) * 255 )
                    #c2 = int( (1.0 - self.fs.curl[self.I(i, j)]) * 255 )
                    if c < 0: c = 0                    
                    color.setRgb(c,c,c)
                    #color.setRgb(c,c2,(c2+c)/2)
                    paint.setPen(color)
                    paint.setBrush(color)
                    paint.drawRect(dx, dy, self.rectSize, self.rectSize)
                    if self.saveImage: self.image.setPixel(i,j,color.rgb())
                # velocity
                if self.drawVelocityField:
                    u = self.fs.u[self.I(i,j)]
                    v = self.fs.v[self.I(i,j)]
                    color.setRgb(255,0,0)
                    paint.setPen(color)
                    dx_mid = dx+rectSize/2
                    dy_mid = dy+rectSize/2
                    paint.drawLine(dx_mid,dy_mid,dx_mid+(u*rectSize),dy_mid+(v*rectSize))
        #self.painter.drawImage(QtCore.QPoint(0,0),self.image) # we could also draw the QImage like this
        paint.end()
        if self.saveImage: self.writeImage()
    
    def paintEvent_image(self, event):
        """ Draw the solver output by filling a QImage, faster than the rect method
            No velocity drawing here for now
            """
        paint = QtGui.QPainter(self)
        # draw grid
        rectSizeDimension = self.rectSize * self.n
        color = QtGui.QColor(255, 255, 255)
        #paint.setPen(color)
        #paint.setBrush(color)
        #paint.drawRect(self.rectSize/2, self.rectSize/2, (self.n-1)*self.rectSize, (self.n-1)*self.rectSize)
        self.image.fill(color.rgb())
        for i in xrange(self.n):
            for j in xrange(self.n):
                dx = (i * (self.rectSize))
                dy = (j * (self.rectSize))
                # density
                if i == 0 or i == self.n-1 or j == 0 or j == self.n-1: continue
                elif self.fs.d[self.I(i, j)] > 0.0085: # making this a small value rather than 0 to avoid drawing those cells
                    # get denisty value from solver density array
                    c = int( (1.0 - self.fs.d[self.I(i, j)]) * 255 )
                    #c2 = int( (1.0 - self.fs.curl[self.I(i, j)]) * 255 )
                    if c < 0: c = 0                    
                    color.setRgb(c,c,c)
                    #color.setRgb(c,c2,(c2+c)/2)
                    #paint.setPen(color)
                    #paint.setBrush(color)
                    #paint.drawRect(dx, dy, self.rectSize, self.rectSize)
                    self.image.setPixel(i,j,color.rgb())
                # velocity
                if self.drawVelocityField:
                    u = self.fs.u[self.I(i,j)]
                    v = self.fs.v[self.I(i,j)]
                    color.setRgb(255,0,0)
                    paint.setPen(color)
                    dx_mid = dx+rectSize/2
                    dy_mid = dy+rectSize/2
                    paint.drawLine(dx_mid,dy_mid,dx_mid+(u*rectSize),dy_mid+(v*rectSize))
        #self.painter.drawImage(QtCore.QPoint(0,0),self.image) # we could also draw the QImage like this
        #paint.end()
        cmpPos = QtCore.QPoint(0,0)
        target = QtCore.QRectF( 0, 0, rectSizeDimension, rectSizeDimension )
        source = QtCore.QRectF( 0, 0, self.n, self.n )
        paint.drawImage(target, self.image, source )
    
    def paintEvent(self, event):
        """ Switching paint event to rects if velocity is requested, should unify that at some point.
            """
        if self.drawVelocityField:
            self.paintEvent_rects(event)
        else:
            self.paintEvent_image(event)    
    
    def writeImage(self):
        self.image.save(self.imagePath+"//"+self.imageName+"."+str(self.frame)+"."+self.imageType.lower(),self.imageType)
        self.frame+=1
    
    def mouseMoveEvent(self, event):
        self.xOld = self.x
        self.yOld = self.y
        self.x = event.pos().x()
        self.y = event.pos().y()
        self.add(event,25)
    
    def mousePressEvent(self, event):
        self.xOld = self.x
        self.yOld = self.y
        self.x = event.pos().x()
        self.y = event.pos().y()
        if event.button() == QtCore.Qt.LeftButton:
            self.add_type = "density"
            self.add(event,75)
        if event.button() == QtCore.Qt.RightButton:
            self.add_type = "velocity"
            self.add(event,75)
    
    def add(self, event, strength=50):
        """ Add sources to solver """
        self.i = int(self.x / self.rectSize)
        self.j = int(self.y / self.rectSize)
        if self.i > self.n: self.i = self.n
        if self.i < 1: self.i = 1
        if self.j > self.n: self.j = self.n
        if self.j < 1: self.j = 1
        if self.add_type == "density":
            # varying the inject size a bit according to the grid size
            injectSize = 1
            if self.n > 69: injectSize = 2
            if self.n > 149: injectSize = 3
            if self.n > 199: injectSize = 4
            for i in xrange(injectSize*-1,injectSize):
                self.fs.dOld[self.I(self.i, self.j+i)] = strength
            for i in xrange(injectSize*-1,injectSize):
                self.fs.dOld[self.I(self.i+i, self.j)] = strength
        if self.add_type == "velocity":
            self.fs.uOld[self.I(self.i,self.j)] = (self.x-self.xOld)/2
            self.fs.vOld[self.I(self.i,self.j)] = (self.y-self.yOld)/2
    
    def keyPressEvent(self, event):
        if type(event) == QtGui.QKeyEvent:
            if event.text() == "r":
                self.fs.reset()
            elif event.text() == "v":
                if not self.drawVelocityField: self.drawVelocityField = True
                else: self.drawVelocityField = False
            elif event.text() == "o":
                if not self.fs.vorticityConfinement_flag: self.fs.vorticityConfinement_flag = True
                else: self.fs.vorticityConfinement_flag = False
            elif event.text() == "b":
                if not self.fs.buoyancy_flag: self.fs.buoyancy_flag = True
                else: self.fs.buoyancy_flag = False
            elif event.text() == ",":
                self.dt -= 0.05
                if self.dt > 1: self.dt = 1
                if self.dt < 0: self.dt = 0
                self.fs.dt = self.dt
            elif event.text() == ".":
                self.dt += 0.05
                if self.dt > 1: self.dt = 1
                if self.dt < 0: self.dt = 0
                self.fs.dt = self.dt
            elif event.text() == "[":
                if not self.n <= 0: self.n-=10
                self.rectSize = self.dimension / self.n
                self.fs.reset(n=self.n)
                self.resize(self.dimension,self.dimension)
                self.image = QtGui.QImage(self.n,self.n,QtGui.QImage.Format_RGB32)
            elif event.text() == "]":
                if not self.n >= 200:self.n+=10
                self.rectSize = self.dimension / self.n
                self.fs.reset(n=self.n)
                self.resize(self.dimension,self.dimension)
                self.image = QtGui.QImage(self.n,self.n,QtGui.QImage.Format_RGB32)
            elif event.text() == "a":
                self.dimension += 20
                self.rectSize = self.dimension / self.n
                self.resize(self.dimension, self.dimension)
            elif event.text() == "z":
                self.dimension -= 20
                self.rectSize = self.dimension / self.n
                self.resize(self.dimension, self.dimension)
            elif event.text() == "d":
                self.fs.diff += 0.00005
            elif event.text() == "c":
                if not self.fs.diff <= 0: self.fs.diff -= 0.00005
                else: self.fs.diff = 0
            elif event.text() == "x":
                if not self.fs.visc <= 0: self.fs.visc -= 0.0005
            elif event.text() == "s":
                self.fs.visc += 0.0005
            elif event.text() == "+" or event.text() == "=":
                if not self.fs.linearSolverIterations>=100: self.fs.linearSolverIterations+=1
            elif event.text() == "-":
                if not self.fs.linearSolverIterations<=1: self.fs.linearSolverIterations-=1
            elif event.text() == "i":
                if self.saveImage: self.saveImage = False
                else: self.saveImage = True
            else:
                event.ignore()
            self.updateStatusBar()
    
    def updateStatusBar(self):
        statusBarText = 'n='+str(self.n)+"x"+str(self.n) + '([ ])  '
        statusBarText += 'dt='+str(self.dt) + '(< >)  '
        statusBarText += 'vortConf='+str(int(self.fs.vorticityConfinement_flag)) + '(o)  '
        statusBarText += 'buoy='+str(int(self.fs.buoyancy_flag)) + '(b)  '
        statusBarText += 'vFld='+str(int(self.drawVelocityField)) + '(v)  '
        statusBarText += 'visc='+str(self.fs.visc) + '(sx)  '
        statusBarText += 'diff='+str(self.fs.diff) + '(dc)  '
        statusBarText += 'linSolvIters='+str(self.fs.linearSolverIterations) + '(+-)  '
        statusBarText += 'writeImage='+str(int(self.saveImage)) + '(i)  '
        self.statusBar().showMessage(statusBarText)
    
    def resizeEvent(self,event):
        self.dimension = self.height()
        self.rectSize = self.dimension / self.n        



class SolverThread(QtCore.QThread):
    def __init__(self, solver):
        QtCore.QThread.__init__(self)
        self.solver = solver
        self.timer = QtCore.QBasicTimer()
        self.timer.start(1/24*1000,self)
    
    def timerEvent(self,event):
        self.solver.velocitySolver()
        self.solver.densitySolver()
        self.emit(QtCore.SIGNAL('solverUpdate()'))   



class FluidSolver(object):
    def __init__(self):
        print "### Python FluidSolver 1.0 - Grid based method as described by Jos Stam. ###"
        self.visc = 0.0
        self.diff = 0.0
        self.linearSolverIterations = 20
        self.vorticityConfinement_flag = True
        self.buoyancy_flag = True
        self.n = 0
        self.dt = 0
        self.size = 0
        self.tmp = []
        self.d = []
        self.dOld = []
        self.u = []
        self.uOld = []
        self.v = []
        self.vOld = []
        self.curl = []
    
    def setup(self, n, dt, viscosity=0.0, diffusion=0.0, vorticityConfinement=True, linearSolverIterations=20):
        self.n = n
        self.dt = dt
        self.size = (self.n+2) * (self.n+2)
        self.visc = viscosity
        self.diff = diffusion
        self.vorticityConfinement_flag = vorticityConfinement
        if linearSolverIterations < 1: linearSolverIterations = 1
        self.linearSolverIterations = linearSolverIterations
        
        self.reset()
    
    def reset(self,n=None):
        """ Reset data structures. """
        if n:
            self.n = n
            self.size = (self.n+2) * (self.n+2)
        self.tmp = [ 0.0 for i in xrange(0,self.size) ]
        self.d = [ 0.0 for i in xrange(0,self.size) ]
        self.dOld = [ 0.0 for i in xrange(0,self.size) ]
        self.u = [ 0.0 for i in xrange(0,self.size) ]
        self.uOld = [ 0.0 for i in xrange(0,self.size) ]
        self.v = [ 0.0 for i in xrange(0,self.size) ]
        self.vOld = [ 0.0 for i in xrange(0,self.size) ]
        self.curl = [ 0.0 for i in xrange(0,self.size) ]
    
    def I(self,i,j):
        return i + (self.n + 2) * j
    
    def buoyancy(self,Fbuoy):
        """ Calculate buoyancy force """
        Tamb = 0
        a = 0.000625
        b = 0.025
        # sum all temperatures
        for i in xrange(1,self.n+1):
            for j in xrange(1,self.n+1):
                Tamb += self.d[self.I(i,j)]   
        # average temperature
        Tamb = Tamb/(self.n*self.n)
        # buoyancy force per cell
        for i in xrange(1,self.n+1):
            for j in xrange(1,self.n+1):
                Fbuoy[self.I(i,j)] = a * self.d[self.I(i,j)] + -b * (self.d[self.I(i,j)] - Tamb)
    
    def calc_curl(self,i,j):
        """ Calculate Curl for cell i, j """
        du_dy = (self.u[self.I(i, j + 1)] - self.u[self.I(i, j - 1)]) * 0.5
        dv_dx = (self.v[self.I(i + 1, j)] - self.v[self.I(i - 1, j)]) * 0.5
        
        return du_dy - dv_dx
    
    def vorticityConfinement(self,Fvc_x,Fvc_y):
        """ Calculate vorticity confinement """
        dw_dx = 0.0
        dw_dy = 0.0
        v = 0.0
        # calculate magnitude of curl(u,v) for each cell
        fabs = math.fabs
        for i in xrange(1,self.n+1):
            for j in xrange(1,self.n+1):
                self.curl[self.I(i,j)] = fabs( self.calc_curl(i, j) )       
        for i in xrange(2, self.n):
            for j in xrange(2, self.n):
                # find derivative of magnitude ( n = del |w| )
                dw_dx = (self.curl[self.I(i + 1, j)] - self.curl[self.I(i - 1, j)]) * 0.5
                dw_dy = (self.curl[self.I(i, j + 1)] - self.curl[self.I(i, j - 1)]) * 0.5
                # calculate vector length (|n|)
                length = float( math.sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001 )
                # N = ( n/|n| )
                dw_dx = dw_dx / length
                dw_dy = dw_dy / length
                # v curl
                v = self.calc_curl(i, j)
                # N x w
                Fvc_x[self.I(i,j)] = dw_dy * -v
                Fvc_y[self.I(i,j)] = dw_dx * v
    
    def velocitySolver(self):
        """ Velocity solving as described by Stam. """
        # add velocity input
        self.addSource(self.u, self.uOld)
        self.addSource(self.v, self.vOld)
        
        # add in vorticity confinement force
        if self.vorticityConfinement_flag:
            self.vorticityConfinement(self.uOld, self.vOld)
            self.addSource(self.u, self.uOld)
            self.addSource(self.v, self.vOld)
        
        # add in buoyancy force
        if self.buoyancy_flag:
            self.buoyancy(self.vOld)
            self.addSource(self.v, self.vOld)
        
        # calculate diffusion in velocity
        self.u, self.uOld = self.uOld, self.u
        self.diffuse(0, self.u, self.uOld, self.visc)
        
        self.v, self.vOld = self.vOld, self.v
        self.diffuse(0, self.v, self.vOld, self.visc)
        
        # create an incompressible field for more effective advection
        self.project(self.u, self.v, self.uOld, self.vOld)
        
        self.u, self.uOld = self.uOld, self.u
        self.v, self.vOld = self.vOld, self.v
        
        # self advect velocities
        self.advect(1, self.u, self.uOld, self.uOld, self.vOld)
        self.advect(2, self.v, self.vOld, self.uOld, self.vOld)
        
        # make an incompressible field
        self.project(self.u, self.v, self.uOld, self.vOld)
        
        # clear all input velocities for next frame
        for i in xrange(self.size):
            self.uOld[i] = 0
            self.vOld[i] = 0
    
    def densitySolver(self):
        """ Density solving """
        # add density input
        self.addSource(self.d, self.dOld)
        self.d, self.dOld = self.dOld, self.d #self.swapD()
        
        self.diffuse(0, self.d, self.dOld, self.diff)
        self.d, self.dOld = self.dOld, self.d #self.swapD()
        
        self.advect(0, self.d, self.dOld, self.u, self.v)
        
        # clear input density array for next frame
        self.dOld = [0 for i in xrange(self.size)]
    
    def addSource(self,x,x0):
        """ Adding sources into data structures """
        for i in xrange(self.size):
            x[i] += self.dt * x0[i]
    
    def advect(self,b,d,d0,du,dv):
        """ Get new velocity for cell, done by backtracing a particle through the velocity field and interpolate the cells were it came from. Stam's method """
        i0 = 0; j0 = 0; i1 = 0; j1 = 0
        x = 0.0; y = 0.0; s0 = 0.0; t0 = 0.0; s1 = 0.0; t1 = 0.0; #dt0 = 0.0
        
        dt0 = self.dt * self.n
        
        for i in xrange(1, self.n+1):
            for j in xrange(1, self.n+1):
                # go backwards through velocity field
                x = i - dt0 * du[self.I(i, j)]
                y = j - dt0 * dv[self.I(i, j)]
                
                # interpolate result
                if x > self.n + 0.5:
                    x = self.n + 0.5
                if x < 0.5:
                    x = 0.5
                
                i0 = int(x)
                i1 = i0 + 1
                
                if y > self.n + 0.5:
                    y = self.n + 0.5
                if y < 0.5:
                    y = 0.5
                
                j0 = int(y)
                j1 = j0 + 1
                
                s1 = x - i0
                s0 = 1 - s1
                t1 = y - j0
                t0 = 1 - t1
                
                d[self.I(i, j)] = s0 * (t0 * d0[self.I(i0, j0)] + t1 * d0[self.I(i0, j1)]) + s1 * (t0 * d0[self.I(i1, j0)] + t1 * d0[self.I(i1, j1)])
        
        self.setBoundry(b, d)
    
    def diffuse(self,b,c,c0,diff):
        a = self.dt * diff * self.n * self.n
        self.linearSolver(b, c, c0, a, 1 + 4 * a)
    
    def project(self, x, y, p, div):
        """ Make the velocity mass conserving, incompressible """
        for i in xrange(1,self.n+1):
            for j in xrange(1,self.n+1):
                div[self.I(i, j)] = (x[self.I(i+1, j)] - x[self.I(i-1, j)] + y[self.I(i, j+1)] - y[self.I(i, j-1)]) * - 0.5 / self.n
                p[self.I(i, j)] = 0
        
        self.setBoundry(0, div)
        self.setBoundry(0, p)
        
        self.linearSolver(0, p, div, 1, 4)
        
        for i in xrange(1,self.n+1):
            for j in xrange(1,self.n+1):
                x[self.I(i, j)] -= 0.5 * self.n * (p[self.I(i+1, j)] - p[self.I(i-1, j)])
                y[self.I(i, j)] -= 0.5 * self.n * (p[self.I(i, j+1)] - p[self.I(i, j-1)])
        
        self.setBoundry(1, x)
        self.setBoundry(2, y)
    
    def linearSolver(self,b,x,x0,a,c):
        """ Iterative linear system solver using the Gauss-sidel relaxation technique. """
        for k in xrange(self.linearSolverIterations):
            for i in xrange(1,self.n+1):
                for j in xrange(1,self.n+1):
                    x[self.I(i, j)] = (a * ( x[self.I(i-1, j)] + x[self.I(i+1, j)] + x[self.I(i, j-1)] + x[self.I(i, j+1)]) + x0[self.I(i, j)]) / c
            self.setBoundry(b, x)
    
    def setBoundry(self,b,x):
        """ set boundary conditions """
        for i in xrange(1,self.n+1):
            if b == 1: x[ self.I(0, i) ] = -x[self.I(1, i)]
            else: x[ self.I(0, i) ] = x[self.I(1, i)]
            if b == 1: x[self.I(self.n+1, i)] = -x[self.I(self.n, i)]
            else: x[self.I(self.n+1, i)] = x[self.I(self.n, i)]
            if b == 2: x[self.I(  i, 0  )] = -x[self.I(i, 1)]
            else: x[self.I(  i, 0  )] = x[self.I(i, 1)]
            if b == 2: x[self.I(  i, self.n+1)] = -x[self.I(i, self.n)]
            else: x[self.I( i, self.n+1)] = x[self.I(i, self.n)]
        
        x[self.I(  0,   0)] = 0.5 * (x[self.I(1, 0  )] + x[self.I(  0, 1)])
        x[self.I(  0, self.n+1)] = 0.5 * (x[self.I(1, self.n+1)] + x[self.I(  0, self.n)])
        x[self.I(self.n+1,   0)] = 0.5 * (x[self.I(self.n, 0  )] + x[self.I(self.n+1, 1)])    
        x[self.I(self.n+1, self.n+1)] = 0.5 * (x[self.I(self.n, self.n+1)] + x[self.I(self.n+1, self.n)])   



class FluidSolverC(FluidSolver):
    """ extending solver with inlined C code called upon via the scipy.weave.inline functionality to enhance solver performance
        """
    def __init__(self):
        FluidSolver.__init__(self)
        print "### C extentions in use ###"
        # I indexer and curl as support code for C inliners
        self.support_code = """
            #include <math.h>
            #include <Python.h>
            int I(int i, int j, int n){
            return (i + (n + 2) * j);
            }
            float calc_curl(int i, int j, int n, PyObject *u, PyObject *v){
            float du_dy = PyFloat_AsDouble(PyList_GetItem(u, I(i, j + 1, n))) - PyFloat_AsDouble(PyList_GetItem(u, I(i, j - 1, n))) * 0.5;
            float dv_dx = PyFloat_AsDouble(PyList_GetItem(v, I(i + 1, j, n))) - PyFloat_AsDouble(PyList_GetItem(v, I(i - 1, j, n))) * 0.5;
            return (float)du_dy - (float)dv_dx;
            }
            """
    
    def reset(self,n=None):
        """ Reset data structures. """
        # todo: convert to numpy arrays (also in C)
        if n:
            self.n = n
            self.size = (self.n+2) * (self.n+2)
        self.tmp = [ 0.0 for i in xrange(0,self.size) ]
        self.d = [ 0.0 for i in xrange(0,self.size) ]
        self.dOld = [ 0.0 for i in xrange(0,self.size) ]
        self.u = [ 0.0 for i in xrange(0,self.size) ]
        self.uOld = [ 0.0 for i in xrange(0,self.size) ]
        self.v = [ 0.0 for i in xrange(0,self.size) ]
        self.vOld = [ 0.0 for i in xrange(0,self.size) ]
        self.curl = [ 0.0 for i in xrange(0,self.size) ]
    #npy.zeros((self.size),npy.float32)#
    
    def buoyancy(self,Fbuoy):
        """ Calculate buoyancy force """
        n = self.n
        d = self.d
        code = """
            float Tamb = 0;
            float a = 0.000625;
            float b = 0.025;
            
            // sum all temperatures
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            Tamb = Tamb + (float)d[I(i, j, n)];
            }
            }
            
            // get average temperature
            Tamb /= (n * n);
            
            // for each cell compute buoyancy force
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            Fbuoy[I(i, j, n)] = a * (float)d[I(i, j, n)] + -b * ((float)d[I(i, j, n)] - Tamb);
            }
            }
            """
        args = ['n','d','Fbuoy']
        weave.inline(code,args,support_code=self.support_code)
    
    def vorticityConfinement(self,Fvc_x,Fvc_y):
        """ Calculate vorticity confinement """
        n = self.n
        u = self.u
        v = self.v
        curl = self.curl
        code = """
            float dw_dx, dw_dy;
            float length;
            float v_curl;
            
            // Calculate magnitude of curl(u,v) for each cell. (|w|)
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            curl[I(i, j, n)] = abs(calc_curl(i, j, n, u, v));
            }
            }
            
            for (int i = 2; i < n; i++)
            {
            for (int j = 2; j < n; j++)
            {
            
            // Find derivative of the magnitude (n = del |w|)
            dw_dx = ((float)curl[I(i + 1, j, n)] - (float)curl[I(i - 1, j, n)]) * 0.5;
            dw_dy = ((float)curl[I(i, j + 1, n)] - (float)curl[I(i, j - 1, n)]) * 0.5;
            
            // Calculate vector length. (|n|)
            // Add small factor to prevent divide by zeros.
            length = (float)sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001;
            
            // N = ( n/|n| )
            dw_dx = dw_dx / length;
            dw_dy = dw_dy / length;
            
            v_curl = calc_curl(i, j, n, u, v);
            
            // N x w
            Fvc_x[I(i, j, n)] = dw_dy * -v_curl;
            Fvc_y[I(i, j, n)] = dw_dx * v_curl;
            }
            }
            """
        args = ['n','u','v','curl','Fvc_x','Fvc_y']
        weave.inline(code,args,support_code=self.support_code)
    
    def addSource(self,x,x0):
        """ Adding sources into data structures """
        size = self.size
        dt = self.dt
        code = """
            for (int i = 0; i < size; i++)
            {
            x[i] = (float)x[i] + ( (float)dt * (float)x0[i] );
            }
            """
        args = ['size','dt','x','x0']
        weave.inline(code,args,support_code=self.support_code)
    
    def advect(self,b,d,d0,du,dv):
        """ Get new velocity for cell, done by backtracing a particle through the velocity field and interpolate the cells were it came from. Stam's method """    
        n = self.n
        dt = self.dt
        code = """
            int i0, j0, i1, j1;
            float x, y, s0, t0, s1, t1, dt0;
            dt0 = dt * n;
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            // go backwards through velocity field
            x = i - dt0 * (float)du[I(i, j,n)];
            y = j - dt0 * (float)dv[I(i, j,n)];
            
            // interpolate results
            if (x > n + 0.5) x = n + 0.5;
            if (x < 0.5)     x = 0.5;
            
            i0 = (int) x;
            i1 = i0 + 1;
            
            if (y > n + 0.5) y = n + 0.5;
            if (y < 0.5)     y = 0.5;
            
            j0 = (int) y;
            j1 = j0 + 1;
            
            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
            
            d[I(i, j,n)] = s0 * (t0 * (float)d0[I(i0, j0,n)] + t1 * (float)d0[I(i0, j1,n)])
            + s1 * (t0 * (float)d0[I(i1, j0,n)] + t1 * (float)d0[I(i1, j1,n)]);
            
            }
            }
            """
        args = ['n','dt','b','d','d0','du','dv']
        weave.inline(code,args,support_code=self.support_code)
        
        self.setBoundry(b, d)
    
    def project(self, x, y, p, div):
        """ Make the velocity a mass conserving, incompressible """
        n = self.n
        code = """
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            div[I(i, j,n)] = ((float)x[I(i+1, j,n)] - (float)x[I(i-1, j,n)]
            + (float)y[I(i, j+1,n)] - (float)y[I(i, j-1,n)])
            * (float)-0.5 / n;
            p[I(i, j,n)] = 0.0;
            }
            }
            """
        args = ['n','x','y','p','div']
        weave.inline(code,args,support_code=self.support_code)
        
        self.setBoundry(0, div)
        self.setBoundry(0, p)
        
        self.linearSolver(0, p, div, 1, 4)
        
        code2 = """
            for (int i = 1; i <= n; i++)
            {
            for (int j = 1; j <= n; j++)
            {
            x[I(i, j,n)] = (float)x[I(i, j,n)] - (0.5 * n * ((float)p[I(i+1, j,n)] - (float)p[I(i-1, j,n)]));
            y[I(i, j,n)] = (float)y[I(i, j,n)] - (0.5 * n * ((float)p[I(i, j+1,n)] - (float)p[I(i, j-1,n)]));
            }
            }
            """
        weave.inline(code2,args,support_code=self.support_code)
        
        self.setBoundry(1, x)
        self.setBoundry(2, y)
    
    def linearSolver(self,b,x,x0,a,c):
        """ Iterative linear system solver using the Gauss-sidel relaxation technique. """
        n = self.n
        # keeping the k-loop in python for setBoundry convenience
        for k in xrange(self.linearSolverIterations):
            code = """
                for (int i = 1; i <= n; i++){
                for (int j = 1; j <= n; j++){
                x[I(i, j,n)] = (a * ( (float)(x[I(i-1, j,n)]) + (float)(x[I(i+1, j,n)]) + (float)(x[I(i, j-1,n)]) + (float)(x[I(i, j+1,n)])) + (float)(x0[I(i, j,n)])) / c;
                }
                }
                """
            args = ['n','b','x','x0','a','c']
            weave.inline(code,args,support_code=self.support_code)
            self.setBoundry(b, x)
    
    def setBoundry(self,b,x):
        """ set boundary conditions """
        n = self.n
        code = """
            for (int i = 1; i <= n; i++)
            {
            x[I(  0, i,n  )] = b == 1 ? -(float)x[I(1, i,n)] : (float)x[I(1, i,n)];
            x[I(n+1, i,n  )] = b == 1 ? -(float)x[I(n, i,n)] : (float)x[I(n, i,n)];
            x[I(  i, 0,n  )] = b == 2 ? -(float)x[I(i, 1,n)] : (float)x[I(i, 1,n)];
            x[I(  i, n+1,n)] = b == 2 ? -(float)x[I(i, n,n)] : (float)x[I(i, n,n)];
            }
            x[I(  0,   0,n)] = 0.5f * ((float)x[I(1, 0,n  )] + (float)x[I(  0, 1,n)]);
            x[I(  0, n+1,n)] = 0.5f * ((float)x[I(1, n+1,n)] + (float)x[I(  0, n,n)]);
            x[I(n+1,   0,n)] = 0.5f * ((float)x[I(n, 0,n  )] + (float)x[I(n+1, 1,n)]);
            x[I(n+1, n+1,n)] = 0.5f * ((float)x[I(n, n+1,n)] + (float)x[I(n+1, n,n)]);
            """
        args = ['n','b','x']
        weave.inline(code,args,support_code=self.support_code)



if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = DrawFluidQt()
    window.show()
    app.exec_()

