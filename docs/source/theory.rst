Mapped-averaging theory
########################

**Absolute free energy**

The Helmholtz free energy :math:`A` of a system at temperature :math:`T` and volume :math:`V` is related to its configurational partition function :math:`Q` via:

.. math::
   A = -k_{\rm B}T \ln{Q}
where :math:`Q` is given by

.. math::
   Q = \int_{V} e^{-\beta U\left({\bf x}\right)} {\rm d} {\bf x}
with :math:`{\bf x}` represents coordinates of all atoms.
For simplicity, we will use unitless energy :math:`{\cal U}\equiv \beta U` and free energy :math:`{\cal A}\equiv \beta A` (and their derivatives; e.g., force).

**Free energy derivatives**

Derivative of free energy w.r.t external distortion (e.g., temperature or volume) is related to material properites. For example, average energy :math:`U` and pressure :math:`P` are given by

.. math::
   U = {\cal A}_{\beta}  \qquad {\rm  and} \qquad  P = -k_{\rm B}T \; {\cal A}_V

To get a general expression of free energy derivative, we will the vector :math:`{\bf \lambda}` to represent all distortions of interest; e.g., :math:`\lambda=\left(\beta, V\right)`. The unitless free energy :math:`{\cal A}` at some :math:`\lambda` is then given by

.. math::
   {\cal A}\left(\lambda\right) = - \ln{Q\left(\lambda\right)}
and its first and second derivatives are given by (with the derivative operators, :math:`\partial_{\nu}\equiv\frac{\partial}{\partial \nu}`)

.. math::
   \partial_{\nu}{\cal A} = -\frac{\partial_{\nu} Q}{Q} \qquad {\rm and} \qquad \partial_{\nu\mu}{\cal A} = -\frac{\partial_{\nu\mu}Q }{Q} + \frac{\partial_{\nu} Q}{Q}  \; \frac{\partial_{\mu} Q}{Q} 

To get :math:`Q` derivatives, we need to recognize that the integration boundary :math:`\Omega` varies with :math:`\lambda`, 

.. math::
   Q\left(\lambda\right) = \int_{\Omega\left(\lambda\right)} e^{-{\cal U}\left({\bf y},\lambda\right)} {\rm d} {\bf y}
where :math:`{\bf y}` represents the new coordinates that span the general phase space :math:`\Omega\left(\lambda\right)`. However, we are only 

.. math::
   = \int_{V} e^{-{\cal U}\left({\bf y}\left({\bf x},\lambda\right),\lambda\right)} J {\rm d} {\bf x}
                         = \int_{V} e^{-{\cal U'}} {\rm d} {\bf x}

where :math:`{\cal U'} \equiv {\cal U} - \ln{J}`

.. math::
   \partial_{\nu} Q &=& - \int_{V}  e^{-{\cal U'}} \; {\rm D}_{\nu} {\cal U'} \;\;  {\rm d}{\bf x}\\
   \partial_{\nu\mu}Q &=& \int_{V} e^{-{\cal U'}}\left[ \left({\rm D}_{\nu} {\cal U'}\right) \left({\rm D}_{\mu} {\cal U'}\right) - {\rm D}_{\nu\mu} {\cal U'} \right] \;  {\rm d}{\bf x}

The first and second free energy derivatives are then given by:

.. math::
   \partial_{\nu}{\cal A} = \; \left< {\rm D}_{\nu} {\cal U'} \right> \qquad {\rm and} \qquad
   \partial_{\nu\mu}{\cal A} = \; \left< {\rm D}_{\nu\mu} {\cal U'} \right>
   - {\rm Cov}\left({\rm D}_{\nu} {\cal U'} \;,\; {\rm D}_{\mu} {\cal U'} \right) 

where, :math:`{\rm Cov}\left(X,Y\right)\equiv \left<XY\right> - \left<X\right> \left<Y\right>`

.. math::
   {\rm D}_{\nu} {\cal U'} &=& {\rm D}_{\nu} {\cal U} - {\rm D}_{\nu} J \\
   {\rm D}_{\nu\mu} {\cal U'} &=& {\rm D}_{\nu\mu} {\cal U} - {\rm D}_{\nu\mu} J 

.. math::
   {\rm D}_{\nu} f = \partial_{\nu} f + {\bf v}^{\nu} \cdot \nabla f

The energy derivatives are given by:

.. math::
   {\rm D}_{\nu} {\cal U} &=& \partial_{\nu} {\cal U} - {\cal F} \cdot {\bf v}^{\nu}\\
   {\rm D}_{\nu\mu}{\cal U} &=& \partial_{\nu\mu} {\cal U} 
   - \left( \partial_{\mu} {\bf v}^{\nu} + {\bf v}^{\mu}\cdot \nabla {\bf v}^{\nu} \right)\cdot {\cal F} 
   + {\bf v}^{\nu} \cdot {\Phi} \cdot {\bf v}^{\mu}
   - \left({\bf v}^{\nu} \cdot \partial_{\mu} {\cal F} 
   + {\bf v}^{\mu} \cdot \partial_{\nu} {\cal F} \right)

.. \left[\left< {\rm D}_{\nu}{\cal U'} {\rm D}_{\mu}{\cal U'} \right> 
.. - \left< {\rm D}_{\nu} {\cal U'} \right>  \left< {\rm D}_{\mu} {\cal U'} \right> 

|

|

|
.. math::
   \partial_{\nu\mu}{\cal A} = \left< {\rm D}_{\nu\mu} {\cal V} \right>
  




   \partial_{\nu}{\cal A} = \left< \partial_{\nu} {\cal U} - \nabla \cdot {\bf v}^{\nu} - {\cal F}\cdot {\bf v}^{\nu}\right>

