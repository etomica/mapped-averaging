Mapped-averaging theory
########################

* **Absolute free energy**

The Helmholtz free energy :math:`A` of a system at temperature :math:`T` and volume :math:`V` is related to its configurational partition function :math:`Q` via:

.. math::
   A = -k_{\rm B}T \ln{Q}
where :math:`Q` is given by

.. math::
   Q = \int_{V} e^{-\beta U\left({\bf x}\right)} {\rm d} {\bf x}
with :math:`{\bf x}` represents coordinates of all atoms.
For simplicity, from now on we will be using unitless energy :math:`{\cal U}\equiv \beta U` and free energy :math:`{\cal A}\equiv \beta A` (and their derivatives; e.g., force).

* **Free energy derivatives**

Derivative of free energy w.r.t external perturpation or distortion (e.g., temperature or volume) is related to material properites. For example, average energy :math:`U` and pressure :math:`P` are given by

.. math::
   U = {\cal A}_{\beta}  \qquad {\rm  and} \qquad  P = -k_{\rm B}T \; {\cal A}_V

To get general expression of free energy derivative, we will use the vector :math:`{\bf \lambda}` to represent all perturpations of interest; e.g., :math:`\lambda=\left(\beta, V\right)`. The unitless free energy :math:`{\cal A}` at some :math:`\lambda` is then given by

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

where we used the :math:`{\rm D}_{\nu}` operator on some function :math:`f({\bf y}({\bf x},\lambda),\lambda)` to represent the total (or, Lagrangian) derivative (i.e., :math:`{\rm D}_{\nu} f = \partial_{\nu} f + \partial_{\nu} {\bf y} \cdot \nabla f`). Accordingly, the derivatives of :math:`{\cal A}` are simply given by

.. math::
   \partial_{\nu}{\cal A} = \; \left< {\rm D}_{\nu} {\cal U'} \right> \qquad {\rm and} \qquad
   \partial_{\mu\nu}{\cal A} = \; \left< {\rm D}_{\mu\nu} {\cal U'} \right>
   - {\rm Cov}\left({\rm D}_{\nu} {\cal U'} \;,\; {\rm D}_{\mu} {\cal U'} \right) 

where, :math:`{\rm Cov}\left(X,Y\right)\equiv \left<XY\right> - \left<X\right> \left<Y\right>` is the covariance between the stochastic variables :math:`X` and :math:`Y`.

* **Evaluation of energy derivatives**

We are left with evaluating derivatives of :math:`{\cal U'}`, which are related to :math:`U` derivatives via

.. math::
   {\rm D}_{\nu} {\cal U'} = {\rm D}_{\nu} {\cal U} - {\rm D}_{\nu} J 
   \qquad {\rm and} \qquad 
   {\rm D}_{\mu\nu} {\cal U'} = {\rm D}_{\mu\nu} {\cal U} - {\rm D}_{\mu\nu} J 


First, :math:`\cal U` derivatives can be directly evaluated using the relation between the total (Lagrangian) and partial (Eulerian) derivatives: 

.. math::
   {\rm D}_{\nu} {\cal U} = \partial_{\nu} {\cal U} - {\cal F} \cdot {\bf v}^{\nu}

where :math:`{\cal F}\equiv -\nabla {\cal U}=-\beta \nabla U` is forces vector on all atoms and :math:`{\bf v}^{\nu}\equiv \partial_{\nu} {\bf y}` represents the mapping "velocity" of the external perturpation :math:`\nu`. Applying this operator twice, we get the second derivative

.. math::
   {\rm D}_{\mu\nu}{\cal U}  = \partial_{\mu\nu} {\cal U} 
   - \left( \partial_{\mu} {\bf v}^{\nu} + {\bf v}^{\mu}\cdot \nabla {\bf v}^{\nu} \right)\cdot {\cal F} 
   + {\bf v}^{\nu} \cdot {\Phi} \cdot {\bf v}^{\mu}
   - \left({\bf v}^{\nu} \cdot \partial_{\mu} {\cal F} 
   + {\bf v}^{\mu} \cdot \partial_{\nu} {\cal F} \right)

where :math:`{\Phi}\equiv \nabla \nabla {\cal U} = \beta \nabla \nabla {\cal U}\;`  is the force constant matrix.

Now, in order to get :math:`{\rm D}_{\nu}J`, we need to do two steps. First, perform diffirentiation of the phase space volume, using (again) the change of variables technique 

.. math::
   \frac{\rm d}{{\rm d}\nu} {\Omega(\lambda)} =
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} 1\; {\rm d} {\bf y} =
   \frac{\rm d}{{\rm d}\nu} \int_{V} {\rm D}_{\nu}J \; {\rm d} {\bf x}

Second,use the Reynolds transport theorem along with the divergence theorem in our multidimentional space

.. math::
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} f\left({\bf y},\lambda\right){\rm d} {\bf y} = \int_{\Omega(\lambda)}     \left[\partial_{\nu} f + \nabla \cdot \left({\bf v}^{\nu} f\right)\right] {\rm d} {\bf y}

Applying this theorem to our case of interest (i.e., :math:`f=1`; hence, :math:`\partial_{\nu}f=0`), we get

.. math::
   \frac{\rm d}{{\rm d}\nu} \int_{\Omega(\lambda)} 1\; {\rm d} {\bf y} = \int_{\Omega(\lambda)}     \nabla \cdot \left({\bf v}^{\nu} f\right) {\rm d} {\bf y}
   =
   \int_{V} \nabla \cdot \left({\bf v}^{\nu} f\right) J {\rm d} {\bf x}

where we used the change of variables in the last term on the right-hand side. Now, equating both derivatives we directly get and expression for :math:`{\rm D}_{\nu}J`

.. math::
   {\rm D}_{\nu}J = J \nabla \cdot {\bf v}^{\nu} 

Repeating the same process with another derivative w.r.t. :math:`\mu`, we directly get

.. math::
   {\rm D}_{\mu\nu}J = J \left[\nabla \cdot \left(\partial_{\mu}{\bf v}^{\nu}\right) 
   + {\bf v}^{\mu}\cdot \nabla\left(\nabla\cdot{\bf v}^{\nu}\right)\right]

Since we are interested at evaluating the derivatives at :math:`{\bf y}={\bf x}`, then :math:`J=1`; hence
:math:`{\rm D}_{\nu}J = \nabla \cdot {\bf v}^{\nu}` and :math:`{\rm D}_{\mu\nu}J = \nabla \cdot \left(\partial_{\mu}{\bf v}^{\nu}\right)  + {\bf v}^{\mu}\cdot \nabla\left(\nabla\cdot{\bf v}^{\nu}\right)`. 


**Final expressions**

.. math::
   \partial_{\nu}{\cal A} = \left< \partial_{\nu} {\cal U} - \nabla \cdot {\bf v}^{\nu} - {\cal F}\cdot {\bf v}^{\nu}\right>

   \left< {\rm D}_{\nu}{\cal U'} {\rm D}_{\mu}{\cal U'} \right> 
   - \left< {\rm D}_{\nu} {\cal U'} \right>  \left< {\rm D}_{\mu} {\cal U'} \right> 

   ......................


**Mapping**

