.. _theory:

Theory
##########

Absolute free energy
====================
The classical Helmholtz configurational free energy :math:`A` of a system at temperature :math:`T` and volume :math:`V` is related to its configurational partition function :math:`Q` via:

.. math::
   A = -k_{\rm B}T \ln{Q} \qquad {\rm where} \qquad 
   Q = \int_{V} e^{-\beta U\left({\bf x}\right)} {\rm d} {\bf x}

with :math:`{\bf x}` representing coordinates of all atoms.
For simplicity, from now on we will be using unitless energy :math:`{\cal U}\equiv \beta U` and free energy :math:`{\cal A}\equiv \beta A` (and their derivatives; e.g., force), where :math:`\beta = 1/k_{\rm B}T`.

Free energy derivatives
=========================

Derivative of free energy w.r.t external perturbation or distortion (e.g., temperature or volume) is related to material properties. For example, average energy :math:`U` and pressure :math:`P` are given by

.. math::
   U = {\cal A}_{\beta}  \qquad {\rm  and} \qquad  P = -k_{\rm B}T \; {\cal A}_V

To get general expression of free energy derivative, we will use the vector :math:`{\bf \lambda}` to represent all perturbations of interest; e.g., :math:`\lambda=\left(\beta, V\right)`. The unitless free energy :math:`{\cal A}` at some :math:`\lambda` is then given by

.. math::
   {\cal A}\left(\lambda\right) = - \ln{Q\left(\lambda\right)}

and its first and second derivatives are given by (using :math:`\partial_{\nu}` to represent the derivative operator :math:`\equiv\frac{\partial}{\partial \nu}`)

.. math::
   \partial_{\nu}{\cal A} = -\frac{\partial_{\nu} Q}{Q} \qquad {\rm and} \qquad \partial_{\mu\nu}{\cal A} = -\frac{\partial_{\mu\nu}Q }{Q} + \frac{\partial_{\nu} Q}{Q}  \; \frac{\partial_{\mu} Q}{Q} 

To get :math:`Q\left(\lambda\right)` derivatives, we need to recognize that, in general, the integration boundary :math:`\Omega` is a function of :math:`\lambda`, 

.. math::
   Q\left(\lambda\right) = \int_{\Omega\left(\lambda\right)} e^{-{\cal U}\left({\bf y},\lambda\right)} {\rm d} {\bf y}

where :math:`{\bf y}` represents the new coordinates that span the general phase space :math:`\Omega\left(\lambda\right)`, at some general :math:`\lambda`. However, we can use the change of variables technique to carry out the integration over some :math:`\lambda`-independent phase-space boundary (we will use :math:`\Omega(\lambda)=V`) and use the Jacobian determinant :math:`J` to transform from the :math:`\bf x` configurations (which span :math:`V`) to the mapped coordinates :math:`\bf y({\bf x},\lambda)`, using :math:`{\rm d}{\bf y} = J {\rm d}{\bf x}` 

.. math::
  Q\left(\lambda\right) = \int_{V} e^{-{\cal U}\left({\bf y}\left({\bf x},\lambda\right),\lambda\right)} J {\rm d} {\bf x}

now, :math:`Q` can be written as a function of the "normal" (current) coordinates :math:`\bf x` 

.. math::
   Q\left(\lambda\right) = \int_{V} e^{-{\cal U'}} {\rm d} {\bf x}

where we define a modified potential :math:`{\cal U'} \equiv {\cal U} - \ln{J}`, which accounts for the mapping.
Now, the :math:`Q` derivatives can be easily evaluated

.. math::
   \partial_{\nu} Q &=& - \int_{V}  e^{-{\cal U'}} \; {\rm D}_{\nu} {\cal U'} \;\;  {\rm d}{\bf x}\\
   \partial_{\mu\nu}Q &=& \int_{V} e^{-{\cal U'}}\left[ \left({\rm D}_{\nu} {\cal U'}\right) \left({\rm D}_{\mu} {\cal U'}\right) - {\rm D}_{\mu\nu} {\cal U'} \right] \;  {\rm d}{\bf x}

where we used the :math:`{\rm D}_{\nu}` operator on some function :math:`f({\bf y}({\bf x},\lambda),\lambda)` to represent the total (Lagrangian or material) derivative (i.e., :math:`{\rm D}_{\nu} f = \partial_{\nu} f + \left(\partial_{\nu} {\bf y}\right) \cdot \nabla f`).  Accordingly, the derivatives of :math:`{\cal A}` are simply given by

.. math::
   \partial_{\nu}{\cal A} = \; \left< {\rm D}_{\nu} {\cal U'} \right> \qquad {\rm and} \qquad
   \partial_{\mu\nu}{\cal A} = \; \left< {\rm D}_{\mu\nu} {\cal U'} \right>
   - {\rm Cov}\left({\rm D}_{\nu} {\cal U'} \;,\; {\rm D}_{\mu} {\cal U'} \right) 

where :math:`{\rm Cov}\left(X,Y\right)\equiv \left<XY\right> - \left<X\right> \left<Y\right>` is the covariance between the stochastic variables :math:`X` and :math:`Y`.
We are left with evaluating derivatives of :math:`{\cal U'}`, which are related to :math:`U` derivatives via :math:`{\rm D}_{\nu} {\cal U'} = {\rm D}_{\nu} {\cal U} - {\rm D}_{\nu} J \; {\rm and} \; {\rm D}_{\mu\nu} {\cal U'} = {\rm D}_{\mu\nu} {\cal U} - {\rm D}_{\mu\nu} J`.

1. *Evaluation of energy derivatives*

First, :math:`\cal U` derivatives can be directly evaluated using the relation between the total (Lagrangian) and partial (Eulerian) derivatives: 

.. math::
   {\rm D}_{\nu} {\cal U} = \partial_{\nu} {\cal U} - {\cal F} \cdot {\dot {\bf x}}^{\nu}

where :math:`{\cal F}\equiv -\nabla {\cal U}=-\beta \nabla U` is the force vector on all atoms and :math:`{\dot {\bf x}}^{\nu}\equiv \partial_{\nu} {\bf y}` represents the mapping "velocity" of the external perturbation :math:`\nu`. Applying this operator twice, we get the second derivative

.. math::
   {\rm D}_{\mu\nu}{\cal U}  = \partial_{\mu\nu} {\cal U} 
   - \left({\ddot {\bf x}}^{\mu\nu} + {\dot {\bf x}}^{\mu}\cdot \nabla {\dot {\bf x}}^{\nu} \right)\cdot {\cal F} 
   + {\dot {\bf x}}^{\nu} \cdot {\Phi} \cdot {\dot {\bf x}}^{\mu}
   - \left({\dot {\bf x}}^{\nu} \cdot \partial_{\mu} {\cal F} 
   + {\dot {\bf x}}^{\mu} \cdot \partial_{\nu} {\cal F} \right)

where :math:`{\Phi}\equiv \nabla \nabla {\cal U} = \beta \nabla \nabla {\cal U}\;`  is the force constant matrix and :math:`{\ddot {\bf x}}^{\mu\nu} \equiv \partial_{\mu} {\dot {\bf x}}^{\nu}` is the :math:`\mu\nu` "acceleration", or the rate of change of :math:`{\dot {\bf x}}^{\nu}` w.r.t. :math:`\mu` (note that :math:`\;{\ddot {\bf x}}^{\mu\nu} \neq {\ddot {\bf x}}^{\nu\mu}`).

|

2. *Evaluation of Jacobian derivatives*

Now, in order to get :math:`{\rm D}_{\nu}J`, we need to do two steps. First, perform differentiation of the phase space volume, using (again) the change of variables technique 

.. math::
   \frac{\rm d}{{\rm d}\nu} {\Omega(\lambda)} =
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} 1\; {\rm d} {\bf y} =
   \frac{\rm d}{{\rm d}\nu} \int_{V} {\rm D}_{\nu}J \; {\rm d} {\bf x}

Second,use the Reynolds transport theorem along with the divergence theorem in our multidimensional space

.. math::
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} f\left({\bf y},\lambda\right){\rm d} {\bf y} = \int_{\Omega(\lambda)}     \left[\partial_{\nu} f + \nabla \cdot \left({\dot {\bf x}}^{\nu} f\right)\right] {\rm d} {\bf y}

Applying this theorem to our case of interest (i.e., :math:`f=1`; hence, :math:`\partial_{\nu}f=0`), we get

.. math::
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} 1\; {\rm d} {\bf y} = \int_{\Omega(\lambda)}     \nabla \cdot \left({\dot {\bf x}}^{\nu} f\right) {\rm d} {\bf y}
   =
   \int_{V} \nabla \cdot \left({\dot {\bf x}}^{\nu} f\right) J {\rm d} {\bf x}

where we used the change of variables in the last term on the right-hand side. Now, equating both derivatives we directly get and expression for :math:`{\rm D}_{\nu}J`

.. math::
   {\rm D}_{\nu}J = J \nabla \cdot {\dot {\bf x}}^{\nu} 

Repeating the same process with another derivative w.r.t. :math:`\mu`, we directly get

.. math::
   {\rm D}_{\mu\nu}J = J \left[\nabla \cdot \left(\partial_{\mu}{\dot {\bf x}}^{\nu}\right) 
   + {\dot {\bf x}}^{\mu}\cdot \nabla\left(\nabla\cdot{\dot {\bf x}}^{\nu}\right)\right]

Since we are interested at evaluating the derivatives at :math:`{\bf y}={\bf x}`, then :math:`J=1`; hence
:math:`{\rm D}_{\nu}J = \nabla \cdot {\dot {\bf x}}^{\nu}` and :math:`{\rm D}_{\mu\nu}J = \nabla \cdot \left(\partial_{\mu}{\dot {\bf x}}^{\nu}\right)  + {\dot {\bf x}}^{\mu}\cdot \nabla\left(\nabla\cdot{\dot {\bf x}}^{\nu}\right)`. 





Mapping velocity
=================
Since :math:`Q` is only a function of :math:`\lambda`, **average** free energy derivatives do not depend on how :math:`{\bf x}` get mapped into the :math:`{\bf y}` coordinates; or, in other words, they do not depend on the mapping velocity :math:`{\dot {\bf x}}^{\nu}`. However, the **fluctuations** (or uncertainty) in these averages do depend on the mapping. Therefore, for the purposes of molecular simulation measurements we need to choose :math:`{\dot {\bf x}^{\nu}}` that reduces the stochastic uncertainty as much as possible.

To develop such a mapping we need to recognize that free energy derivatives are given as ensemble averages over :math:`{\rm D}_{\nu} {\cal U'}` (and its derivative, :math:`{\rm D}_{\mu\nu} {\cal U'}`).
Therefore, a perfect mapping is such that :math:`{\rm D}_{\nu} {\cal U'}` is independent on coordinates :math:`\bf x`; hence

.. math::
   \partial_{\nu}{\cal A} = \; \left< {\rm D}_{\nu} {\cal U'} \right> 
   = {\rm D}_{\nu} {\cal U'}

Using the above energy and Jacobian derivatives, we get

.. math::
   \partial_{\nu}{\cal A} = \partial_{\nu} {\cal U} - \nabla \cdot {\dot {\bf x}}^{\nu} - {\cal F}\cdot {\dot {\bf x}}^{\nu}

Solving this equation yields the unique mapping that yields no fluctuations; however, there are two problems. First of all, :math:`\partial_{\nu}{\cal A}` is the very quantity that we need to measure. Second, since :math:`{\dot {\bf x}}^{\nu}` is a multidimensional vector (:math:`3N` for the case of atomic systems) we have under-determined system as we only have one equation to solve. 

The first problem is solved using the fast that :math:`{\dot {\bf x}}^{\nu}` does not affect average estimates; hence, it can be derived from another (known) system, which we will call reference. 

.. math::
   \partial_{\nu}{\cal A}^{\rm ref} = \partial_{\nu} {\cal U}^{\rm ref} - \nabla \cdot {\dot {\bf x}}^{\nu} - {\cal F}^{\rm ref}\cdot {\dot {\bf x}}^{\nu}

where :math:`\partial_{\nu}{\cal A}^{\rm ref}` is a reference-dependent constant (only function of :math:`\lambda`), named :math:`c`.

To solve the second problem, we will assume that each degree of freedom (dof) is mapped with the same amount (scaling); so

.. math::
   \partial_{\nu} {\cal u}^{\rm ref} - \partial_{x} {\dot x}^{\nu} - {\cal f}^{\rm ref} {\dot x}^{\nu} = \partial_{\nu}{\cal a}^{\rm ref} \equiv c(\lambda) 

where small symbols represent an intensive quantities (i.e., :math:`x\equiv X/{\rm dof}`). For a given :math:`\lambda`, this is a standard first-order differential equation, with the unknown being the velocity of mapping :math:`{\dot x}(x,\lambda)`. For simplicity, we will drop the :math:`\lambda` dependency from all terms, hence

.. math::
    \partial_{x} {\dot x}^{\nu}\left( x\right) + {\cal f}\left( x\right)^{\rm ref} {\dot x}^{\nu}\left( x\right)  =
    \partial_{\nu}{\cal u}\left( x\right)^{\rm ref} - c \equiv g\left( x\right) 

where :math:`g(x)` is a known function once a reference system is chosen. The solution of this equation is given by

.. math::
   {\dot x}^{\nu} = e^{-I(x)} \left(\int g \; e^{I(x)}{\rm d}x + {\rm constant} \right) 

where :math:`I(x) \equiv \int f(x)^{\rm ref} {\rm d}x`. The integration constant can be evaluated by requiring the mapping to have some value at some coordinate :math:`x`.


