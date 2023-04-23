Download Link: https://assignmentchef.com/product/solved-aero4630-project1-method-of-manufactured-solutions
<br>
<h2>Part 1a: Horizontal Change</h2>

For and f(x) = −6x<sub>1</sub>, the relationship

−∇<sup>2</sup>u(x) = f(x),x ∈ ω                                                         (1)

can be shown to be true by taking the Laplacian of u(x).

The Python script used to discretize and solve the relationship is shown <strong>Appendix 1</strong>. The Paraview output of is shown in Figure 1.

The error reported by the script output is shown below:

errorL2   = 0.004931859322266561 errormax = 3.33066907388e−16

<h2>Part 1b: Vertical Change</h2>

This relationship can be shown for and f(x) = −6x<sub>2 </sub>using the same steps as were used for and f(x) = −6x<sub>1</sub>. The error from the python script used was the same as well. The Paraview output is shown in Figure 2, and the Python script used is given in <strong>Appendix 1b</strong>.

Figure 1: Paraview output for problem 1a

Figure 2: Paraview output for problem 1b

<h1>Problem 2: Improving Accuracy</h1>

<h2>Part 2a: Improving Through Mesh Resolution</h2>

The error was reduced by increasing the mesh resolution in the python script. Mesh resolutions of 4, 8, 16, 32, and 64 produced different errors. These are shown in Figure 3.

Figure 3: Error vs Mesh Sizing for Problem 2a

The code used to assess mesh refinement is provided in <strong>Appendix 2a</strong>.

<h2>Part 2b: Improving Through the Shape Function</h2>

The error was reduced by increasing the order of or the shape functions. Orders of 1, 2, 3, 4, and 5 were used. A plot of the error produced over the shape function order is given in Figure 4

The code used to assess shape function order is provided in <strong>Appendix 2b</strong>.

Figure 4: Error vs Mesh Sizing for Problem 2a

<h1>Appendix 1a: Code for Problem 1</h1>

”””

Python        script      for    1a   of     Project     1

Original         Author :                       Vinamra Agrawal

Date :                                                                    January    25 ,    2019

Edited By:                                                             Omkar Mulekar

Date :                                                                    January    30 ,    2019

−Laplace (u) = f in the unit square u = u D on the boundary

In this problem , u D = xˆ3

f = −6x

”””

#============================================================

# This         statement      needs    to    be    added    for     most ( i f      not    a l l )

# scripts          for       Fenics <strong>from    </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib . pyplot as        plt

#============================================================

#============================================================ # This    is          to create  a mesh .

# In the class we saw how our rectangular domain was s p l i t # into triangles . This function , by defualt , creates a box # of dimensions 1×1 , with 8 intervals on each edge .

# These intervals are connected through triangles . mesh = UnitSquareMesh (8 , 8)

#============================================================

#============================================================ # This    is          to choose what   kind     of         function          space we are   dealing # with .          Essentially       these are       shape   functions .

# Don ’ t       worry ,    we   w i l l     do shape        functions      in     class       later .

V = FunctionSpace (mesh ,               ’P’ , 1)

#============================================================

<em> </em>

# This         statement       defines      the     boundary       condition     u D

# ”Expression” is a convenient way to declare complicated # expressions in Fenics . On the back−end , a program reads this # and converts i t in a format , the code can work with .

# This is a feature to make your l i f e easier . # degree = 2 t e l l s Fenics , the degree of the polynomial .

u D = Expression ( ’x [0]∗ x [0]∗ x [0] ’ , degree=3)

#============================================================

#============================================================

# This        is       a ” function ”        called      boundary ,     that     takes     in   two     inputs

#                                    x                                                         − location       of   a    point

#                                         on boundary − location         of    the     boundary

# Right now,            this       function      is    not      doing much.      It    is     simply

# returning the location of the boundary . In the future , we w i l l # mess with this function as well . We w i l l see that i t is # extremely useful .

<strong>def </strong>boundary(x ,           on boundary ):

<strong>return </strong>on boundary

#============================================================

#============================================================

# This        is     where we      actually       define     the    boundary       condition

# Fenics only needs you to define Dirichlet ( or displacement ) # boundary condition . You define the Neumann ( or force ) # boundary condition in a different way .

# This         function      takes     in     three      inputs

#                 V                                                       − the       function      space     ( shape         functions )

#                 u D                               − the       expression      of    the   d i r i c h l e t    boundary

#                  boundary                     − a      function       ( declared      above )    that    t e l l s   you

#         exactly where the boundary is bc = DirichletBC (V, u D , boundary)

#============================================================

#============================================================

<table width="590">

 <tbody>

  <tr>

   <td width="191"># This        is          where you</td>

   <td width="65">define</td>

   <td width="196">the      weak form     of       the</td>

   <td width="137">equation</td>

  </tr>

  <tr>

   <td width="191"># First          step :          t e l l</td>

   <td width="65">Fenics</td>

   <td width="196">to       start            with u as a</td>

   <td width="137">t r i a l         function</td>

  </tr>

  <tr>

   <td width="191">#u = TrialFunction (V)</td>

   <td width="65"> </td>

   <td width="196">                in    the          function</td>

   <td width="137">space V</td>

  </tr>

  <tr>

   <td width="191"># Second        step :        t e l l</td>

   <td width="65">Fenics</td>

   <td width="196">to       choose     v(x)         as a</td>

   <td width="137">test        function</td>

  </tr>

  <tr>

   <td width="191">#</td>

   <td width="65"> </td>

   <td width="196">                in    the          function</td>

   <td width="137">space V.</td>

  </tr>

 </tbody>

</table>

v = TestFunction (V)

# Third step : Define f (x) f = Expression ( ’−6∗x [0] ’ , degree=1) # Fourth  step :   Define  the       l e f t    hand           side      of         the       weak form a = dot( grad (u) , grad (v))∗dx

# Fifth step : Define the right hand side of the weak form L = f∗v∗dx

#============================================================

#============================================================ # This    chunk  of code     computes        the       solution .

# First , move away from t r i a l function u and make i t a function . # Basically i t t e l l s Fenics to make u as a combination of # shape functions . We’ l l do this later too . u = Function (V)

# Finally . . . l e t i t run . Solve LHS = RHS, for u given BC. solve (a == L, u, bc)

#============================================================

#============================================================ # Save   solution to         f i l e    in VTK format .

# In       the     class ,     we      visualized       the    f i l e    in     Paraview .

vtkfile = File ( ’ poisson/ solution . pvd ’ ) vtkfile &lt;&lt; u

#============================================================

#============================================================

# Now we       test     i f    our    code    is     working .     This     is       where we use      the

# method       of     manufactured          solutions . We know        that

#                        u = 1 + xˆ2 + 2yˆ2

# s a t i s f i e s the PDE. That ’ s why we chose this as the boundary # condition . u = u D is automatically s a t i s f i e d .

# This        line      computes        the L2 error        between     the    computed      solution

# u and the the solution we expect u D error L2 = errornorm (u D , u, ’L2 ’ )

# Next , we find the value of the solution ( expected and computed) # at each vertex of the mesh . Then we compare the two . Our goal here # is to find the maximum error we get .

vertex values u D = u D. compute vertex values (mesh) vertex values u = u. compute vertex values (mesh)

<strong>import </strong>numpy as np

error max = np.<strong>max</strong>(np. <strong>abs</strong>( vertex values u D − vertex values u ))

# Print errors <strong>print</strong>( ’ error L2 =’ , errorL2 ) <strong>print</strong>( ’ error max =’ , errormax )




<h1>Appendix 1b: Code for Problem 1</h1>

”””

Python        script      for    1b    of     Project     1

Original         Author :                       Vinamra Agrawal

Date :                                                                    January    25 ,    2019

Edited By:                                                             Omkar Mulekar

Date :                                                                    January    30 ,    2019

−Laplace (u) = f in the unit square u = u D on the boundary

In this problem , u D = xˆ3

f = −6x

”””

#============================================================

# This         statement      needs    to    be    added    for     most ( i f      not    a l l )

# scripts          for       Fenics <strong>from    </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib . pyplot as            plt

#============================================================

#============================================================ # This    is          to create  a mesh .

# In the class we saw how our rectangular domain was s p l i t # into triangles . This function , by defualt , creates a box # of dimensions 1×1 , with 8 intervals on each edge .

# These intervals are connected through triangles . mesh = UnitSquareMesh (8 , 8)

#============================================================

#============================================================ # This    is          to choose what   kind     of         function          space we are   dealing # with .          Essentially these   are       shape   functions .

# Don ’ t       worry ,    we   w i l l     do shape        functions      in     class       later .

V = FunctionSpace (mesh ,               ’P’ , 1)

#============================================================ # This    statement defines the       boundary        condition        u D

# ”Expression” is a convenient way to declare complicated # expressions in Fenics . On the back−end , a program reads this # and converts i t in a format , the code can work with .

# This is a feature to make your l i f e easier . # degree = 2 t e l l s Fenics , the degree of the polynomial .

u D = Expression ( ’x [1]∗ x [1]∗ x [1] ’ , degree=3)

#============================================================

#============================================================

# This        is       a ” function ”        called      boundary ,     that     takes     in   two     inputs

#                                    x                                                         − location       of   a    point

#                                         on boundary − location         of    the     boundary

# Right now,            this       function      is    not      doing much.      It    is     simply

# returning the location of the boundary . In the future , we w i l l # mess with this function as well . We w i l l see that i t is # extremely useful .

<strong>def </strong>boundary(x ,           on boundary ):

<strong>return </strong>on boundary

#============================================================

#============================================================

# This        is     where we      actually       define     the    boundary       condition

# Fenics only needs you to define Dirichlet ( or displacement ) # boundary condition . You define the Neumann ( or force ) # boundary condition in a different way .

# This         function      takes     in     three      inputs

#                 V                                                       − the       function      space     ( shape         functions )

#                 u D                               − the       expression      of    the   d i r i c h l e t    boundary

#                  boundary                     − a      function       ( declared      above )    that    t e l l s   you

#         exactly where the boundary is bc = DirichletBC (V, u D , boundary)

#============================================================

#============================================================

<table width="590">

 <tbody>

  <tr>

   <td width="191"># This        is          where you</td>

   <td width="65">define</td>

   <td width="196">the      weak form     of       the</td>

   <td width="137">equation</td>

  </tr>

  <tr>

   <td width="191"># First          step :          t e l l</td>

   <td width="65">Fenics</td>

   <td width="196">to       start            with u as a</td>

   <td width="137">t r i a l         function</td>

  </tr>

  <tr>

   <td width="191">#u = TrialFunction (V)</td>

   <td width="65"> </td>

   <td width="196">                in    the          function</td>

   <td width="137">space V</td>

  </tr>

  <tr>

   <td width="191"># Second        step :        t e l l</td>

   <td width="65">Fenics</td>

   <td width="196">to       choose     v(x)         as a</td>

   <td width="137">test        function</td>

  </tr>

  <tr>

   <td width="191">#</td>

   <td width="65"> </td>

   <td width="196">                in    the          function</td>

   <td width="137">space V.</td>

  </tr>

 </tbody>

</table>

v = TestFunction (V)

# Third step : Define f (x) f = Expression ( ’−6∗x [1] ’ , degree=1) # Fourth  step :   Define  the       l e f t    hand           side      of         the       weak form a = dot( grad (u) , grad (v))∗dx

# Fifth step : Define the right hand side of the weak form L = f∗v∗dx

#============================================================

#============================================================ # This    chunk  of code     computes        the       solution .

# First , move away from t r i a l function u and make i t a function . # Basically i t t e l l s Fenics to make u as a combination of # shape functions . We’ l l do this later too . u = Function (V)

# Finally . . . l e t i t run . Solve LHS = RHS, for u given BC. solve (a == L, u, bc)

#============================================================

#============================================================ # Save   solution to         f i l e    in VTK format .

# In       the     class ,     we      visualized       the    f i l e    in     Paraview .

vtkfile = File ( ’ poisson/ solution . pvd ’ ) vtkfile &lt;&lt; u

#============================================================

#============================================================

# Now we       test     i f    our    code    is     working .     This     is       where we use      the

# method       of     manufactured          solutions . We know        that

#                        u = 1 + xˆ2 + 2yˆ2

# s a t i s f i e s the PDE. That ’ s why we chose this as the boundary # condition . u = u D is automatically s a t i s f i e d .

# This        line      computes        the L2 error        between     the    computed      solution

# u and the the solution we expect u D error L2 = errornorm (u D , u, ’L2 ’ )

# Next , we find the value of the solution ( expected and computed) # at each vertex of the mesh . Then we compare the two . Our goal here # is to find the maximum error we get .

vertex values u D = u D. compute vertex values (mesh) vertex values u = u. compute vertex values (mesh)

<strong>import </strong>numpy as np

error max = np.<strong>max</strong>(np. <strong>abs</strong>( vertex values u D − vertex values u ))

# Print errors <strong>print</strong>( ’ error L2 =’ , errorL2 ) <strong>print</strong>( ’ error max =’ , errormax )




<h1>Appendix 2a: Code for Problem 2a</h1>

<table width="447">

 <tbody>

  <tr>

   <td width="224">”””</td>

   <td width="223"> </td>

  </tr>

  <tr>

   <td width="224">Python       script      for    2a      of</td>

   <td width="223">Project       1</td>

  </tr>

  <tr>

   <td width="224">Original        Author :</td>

   <td width="223">Vinamra Agrawal</td>

  </tr>

  <tr>

   <td width="224">Date :</td>

   <td width="223">                           January    25 ,    2019</td>

  </tr>

  <tr>

   <td width="224">Edited By:</td>

   <td width="223">Omkar Mulekar</td>

  </tr>

  <tr>

   <td width="224">Date :</td>

   <td width="223">                           January    30 ,    2019</td>

  </tr>

  <tr>

   <td width="224">−Laplace (u) = f</td>

   <td width="223">          in    the     unit     square</td>

  </tr>

  <tr>

   <td width="224">                           u = u D       on      the</td>

   <td width="223">boundary</td>

  </tr>

 </tbody>

</table>

In this problem , u D = xˆ3

f = −6x

”””

#============================================================

# This         statement      needs    to    be    added    for     most ( i f      not    a l l )

# scripts          for       Fenics <strong>from    </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use ( ’Agg ’ )

<strong>import </strong>matplotlib . pyplot as                  plt

#============================================================

#============================================================ # This    is          to create  a mesh .

# In the class we saw how our rectangular domain was s p l i t # into triangles . This function , by defualt , creates a box # of dimensions 1×1 , with 8 intervals on each edge . # These intervals are connected through triangles .

# Define Mesh Sizes and i n i t i a l i z e error variables meshSizes = [4 ,8 ,16 ,32 ,64] errorL2 = [0] ∗ <strong>len</strong>( meshSizes ) errormax = [0] ∗ <strong>len</strong>( meshSizes )

<strong>for </strong>i <strong>in range</strong>(<strong>len</strong>( meshSizes )):

mesh = UnitSquareMesh( meshSizes [ i ] ,           meshSizes [ i ])

<em> </em>

#============================================================

#============================================================ # This     is          to         choose what   kind     of         function          space we are   dealing # with .        Essentially       these   are       shape   functions .

# Don ’ t     worry ,    we   w i l l     do shape       functions      in     class       later .

V = FunctionSpace (mesh ,           ’P’ , 1)

#============================================================

#============================================================

# This       statement       defines      the     boundary       condition     u D

# ”Expression” is a convenient way to declare complicated # expressions in Fenics . On the back−end , a program reads this # and converts i t in a format , the code can work with .

# This is a feature to make your l i f e easier . # degree = 2 t e l l s Fenics , the degree of the polynomial .

u D = Expression ( ’x [0]∗ x [0]∗ x [0] ’ , degree=3)

#============================================================

#============================================================

# This      is       a ” function ”        called      boundary ,     that     takes     in   two     inputs

#                                    x                                                         − location       of   a    point

#                                         on boundary − location         of    the     boundary

# Right now,        this       function      is    not      doing much.      It    is     simply

# returning the location of the boundary . In the future , we w i l l # mess with this function as well . We w i l l see that i t is # extremely useful .

<strong>def </strong>boundary(x ,        on boundary ):

<strong>return </strong>on boundary

#============================================================

#============================================================

# This      is     where we      actually       define     the     boundary       condition

# Fenics only needs you to define Dirichlet ( or displacement ) # boundary condition . You define the Neumann ( or force ) # boundary condition in a different way .

# This       function      takes     in     three      inputs

# V − the function space ( shape functions # u D − the expression of the d i r i c h l e t boundary

#                  boundary                     − a      function       ( declared      above )    that    t e l l s   you

#   exactly where the boundary is bc = DirichletBC (V, u D , boundary)

#============================================================

#============================================================

<table width="590">

 <tbody>

  <tr>

   <td width="256"># This        is      where you           define</td>

   <td width="196">the      weak form     of       the</td>

   <td width="137">equation</td>

  </tr>

  <tr>

   <td width="256"># First          step :         t e l l       Fenics</td>

   <td width="196">to       start            with u as a</td>

   <td width="137">t r i a l         function</td>

  </tr>

  <tr>

   <td width="256">#u = TrialFunction (V)</td>

   <td width="196">                in    the          function</td>

   <td width="137">space V</td>

  </tr>

  <tr>

   <td width="256"># Second        step :       t e l l       Fenics</td>

   <td width="196">to       choose     v(x)         as a</td>

   <td width="137">test        function</td>

  </tr>

  <tr>

   <td width="256">#</td>

   <td width="196">                in    the          function</td>

   <td width="137">space V.</td>

  </tr>

 </tbody>

</table>

v = TestFunction (V)

# Third step : Define f (x) f = Expression ( ’−6∗x [0] ’ , degree=1)

# Fourth  step :   Define  the       l e f t    hand    side      of         the       weak form a = dot( grad (u) , grad (v))∗dx

# Fifth step : Define the right hand side of the weak form L = f∗v∗dx

#============================================================

#============================================================ # This     chunk  of         code     computes        the       solution .

# First , move away from t r i a l function u and make i t a function . # Basically i t t e l l s Fenics to make u as a combination of # shape functions . We’ l l do this later too . u = Function (V)

# Finally . . . l e t i t run . Solve LHS = RHS, for u given BC. solve (a == L, u, bc)

#============================================================

#============================================================ # Save     solution           to         f i l e    in VTK format .

# In     the      class ,     we      visualized       the    f i l e    in     Paraview .

vtkfile = File ( ’ poisson/ solution . pvd ’ ) vtkfile &lt;&lt; u

#============================================================

#============================================================

# Now we     test     i f    our    code    is     working .     This     is      where we use      the

# method     of     manufactured          solutions . We know        that

#                        u = 1 + xˆ2 + 2yˆ2

# s a t i s f i e s the PDE. That ’ s why we chose this as the boundary # condition . u = u D is automatically s a t i s f i e d .

# This      line     computes        the L2 error        between     the    computed      solution

# u and the the solution we expect u D error L2 [ i ] = errornorm (u D , u, ’L2 ’ )

# Next ,     we    find     the     value     of    the      solution       ( expected       and computed)

# at           each    vertex  of         the       mesh .  Then we compare       the       two .    Our goal     here # is          to         find      the maximum error we get .

vertex values u D = u D. compute vertex values (mesh) vertex values u = u. compute vertex values (mesh)

<strong>import </strong>numpy as np

error max [ i ] = np.<strong>max</strong>(np. <strong>abs</strong>( vertex values u D − vertex values u ))

# Print errors <strong>print</strong>( ’ meshSizes =’ , meshSizes ) <strong>print</strong>( ’ error L2 =’ , errorL2 ) <strong>print</strong>( ’ error max =’ , errormax )

#============================================================

plt . figure (1) plt . subplot (211) plt . plot ( meshSizes , error L2 , ’bo−’ ) plt . ylabel ( ’ error L2 ’ ) plt . subplot (212) plt . plot ( meshSizes , error max , ’bo−’ ) plt . xlabel ( ’MeshSizes ’ ) plt . ylabel ( ’ errormax ’ ) plt . savefig ( ’2 afig1 . png ’ )

<h1>Appendix 2a: Code for Problem 2b</h1>

<table width="447">

 <tbody>

  <tr>

   <td width="224">”””</td>

   <td width="223"> </td>

  </tr>

  <tr>

   <td width="224">Python       script      for    2b      of</td>

   <td width="223">Project       1</td>

  </tr>

  <tr>

   <td width="224">Original        Author :</td>

   <td width="223">Vinamra Agrawal</td>

  </tr>

  <tr>

   <td width="224">Date :</td>

   <td width="223">                           January    25 ,    2019</td>

  </tr>

  <tr>

   <td width="224">Edited By:</td>

   <td width="223">Omkar Mulekar</td>

  </tr>

  <tr>

   <td width="224">Date :</td>

   <td width="223">                           January    30 ,    2019</td>

  </tr>

  <tr>

   <td width="224">−Laplace (u) = f</td>

   <td width="223">          in    the     unit     square</td>

  </tr>

  <tr>

   <td width="224">                           u = u D       on      the</td>

   <td width="223">boundary</td>

  </tr>

 </tbody>

</table>

In this problem , u D = xˆ3

f = −6x

”””

#============================================================

# This         statement      needs    to    be    added    for     most ( i f      not    a l l )

# scripts          for       Fenics <strong>from    </strong>future            <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use ( ’Agg ’ )

<strong>import </strong>matplotlib . pyplot as                  plt

#============================================================

#============================================================ # This is          to         create  a mesh .

# In the class we saw how our rectangular domain was s p l i t # into triangles . This function , by defualt , creates a box # of dimensions 1×1 , with 8 intervals on each edge . # These intervals are connected through triangles .

# Define Mesh Sizes and i n i t i a l i z e error variables order = [1 ,2 ,3 ,4 ,5] errorL2 = [0] ∗ <strong>len</strong>( order ) errormax = [0] ∗ <strong>len</strong>( order )

<strong>for </strong>i <strong>in range</strong>(<strong>len</strong>( order )):

mesh = UnitSquareMesh (8 , 8) #============================================================

#============================================================ # This     is          to         choose what   kind     of         function          space we are   dealing # with .        Essentially       these   are       shape   functions .

# Don ’ t     worry ,    we   w i l l     do shape       functions      in     class       later .

V = FunctionSpace (mesh ,            ’P’ , order [ i ])

#============================================================

#============================================================

# This       statement       defines      the     boundary       condition     u D

# ”Expression” is a convenient way to declare complicated # expressions in Fenics . On the back−end , a program reads this # and converts i t in a format , the code can work with .

# This is a feature to make your l i f e easier . # degree = 2 t e l l s Fenics , the degree of the polynomial .

u D = Expression ( ’x [0]∗ x [0]∗ x [0] ’ , degree=3)

#============================================================

#============================================================

# This      is       a ” function ”        called      boundary ,     that     takes     in   two     inputs

#                                    x                                                         − location       of   a    point

#                                         on boundary − location         of    the     boundary

# Right now,        this       function      is    not      doing much.      It    is     simply

# returning the location of the boundary . In the future , we w i l l # mess with this function as well . We w i l l see that i t is # extremely useful .

<strong>def </strong>boundary(x ,        on boundary ):

<strong>return </strong>on boundary

#============================================================

#============================================================

# This      is     where we      actually       define     the     boundary       condition

# Fenics only needs you to define Dirichlet ( or displacement ) # boundary condition . You define the Neumann ( or force ) # boundary condition in a different way .

# This       function      takes     in     three      inputs

# V − the function space ( shape functions # u D − the expression of the d i r i c h l e t boundary

#                  boundary                     − a      function       ( declared      above )    that    t e l l s   you

#   exactly where the boundary is bc = DirichletBC (V, u D , boundary)

#============================================================

#============================================================

<table width="590">

 <tbody>

  <tr>

   <td width="256"># This        is      where you           define</td>

   <td width="196">the      weak form     of       the</td>

   <td width="137">equation</td>

  </tr>

  <tr>

   <td width="256"># First          step :         t e l l       Fenics</td>

   <td width="196">to       start            with u as a</td>

   <td width="137">t r i a l         function</td>

  </tr>

  <tr>

   <td width="256">#u = TrialFunction (V)</td>

   <td width="196">                in    the          function</td>

   <td width="137">space V</td>

  </tr>

  <tr>

   <td width="256"># Second        step :       t e l l       Fenics</td>

   <td width="196">to       choose     v(x)         as a</td>

   <td width="137">test        function</td>

  </tr>

  <tr>

   <td width="256">#</td>

   <td width="196">                in    the          function</td>

   <td width="137">space V.</td>

  </tr>

 </tbody>

</table>

v = TestFunction (V)

# Third step : Define f (x) f = Expression ( ’−6∗x [0] ’ , degree=1)

# Fourth  step :   Define  the       l e f t    hand    side      of         the       weak form a = dot( grad (u) , grad (v))∗dx

# Fifth step : Define the right hand side of the weak form L = f∗v∗dx

#============================================================

#============================================================ # This     chunk  of         code     computes        the       solution .

# First , move away from t r i a l function u and make i t a function . # Basically i t t e l l s Fenics to make u as a combination of # shape functions . We’ l l do this later too . u = Function (V)

# Finally . . . l e t i t run . Solve LHS = RHS, for u given BC. solve (a == L, u, bc)

#============================================================

#============================================================ # Save     solution           to         f i l e    in VTK format .

# In     the      class ,     we      visualized       the    f i l e    in     Paraview .

vtkfile = File ( ’ poisson/ solution . pvd ’ ) vtkfile &lt;&lt; u

#============================================================

#============================================================

# Now we     test     i f    our    code    is     working .     This     is      where we use      the

# method     of     manufactured          solutions . We know        that

#                        u = 1 + xˆ2 + 2yˆ2

# s a t i s f i e s the PDE. That ’ s why we chose this as the boundary # condition . u = u D is automatically s a t i s f i e d .

# This      line     computes        the L2 error        between     the    computed      solution

# u and the the solution we expect u D error L2 [ i ] = errornorm (u D , u, ’L2 ’ )

# Next ,     we    find     the     value     of    the      solution       ( expected       and computed)

# at           each    vertex  of         the       mesh .  Then we compare       the       two .    Our goal     here # is          to         find      the maximum error we get .

vertex values u D = u D. compute vertex values (mesh) vertex values u = u. compute vertex values (mesh)

<strong>import </strong>numpy as np

error max [ i ] = np.<strong>max</strong>(np. <strong>abs</strong>( vertex values u D − vertex values u ))

# Print errors <strong>print</strong>( ’ order =’ , order ) <strong>print</strong>( ’ error L2 =’ , errorL2 ) <strong>print</strong>( ’ error max =’ , errormax )

#============================================================

plt . figure (1) plt . subplot (211) plt . plot ( order , error L2 , ’bo−’ ) plt . ylabel ( ’ error L2 ’ ) plt . subplot (212) plt . plot ( order , error max , ’bo−’ ) plt . xlabel ( ’Mesh Sizes ’ ) plt . ylabel ( ’ error max ’ ) plt . savefig ( ’2bfig1 . png ’ )