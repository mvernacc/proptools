Electric Propulsion Basics
**************************

Electrostatic Thrust Generation
===============================
Like all rockets, electric thrusters generate thrust by expelling mass at high velocity. In electric (electrostatic) propulsion, a high-velocity jet is produced by accelerating charged particles through an electric potential.

The propellant medium in electrostatic thrusters is a plasma consisting of electrons and positively charged ions. The plasma is usually produced by ionizing a gas via electron bombardment. The electric fields within the thruster are configured such that there is a large potential drop between the plasma generation region and the thruster exit. When ions exit the thruster, the are accelerated to high velocities by the electrostatic force. If the ions fall through a beam potential of :math:`V_b`, they will leave the thruster with a velocity

.. math::

  v_e = \sqrt{2 V_b \frac{q}{m_i}}

where :math:`q` is the ion charge and  :math:`m_i` is the ion mass. The thrust produced by the ion flow is

.. math::

  F &= v_e \dot{m}

    &= \left( \sqrt{2 V_b \frac{q}{m_i}} \right) \left( \frac{m_i}{q} I_b \right)

    &= I_b \sqrt{2 V_b \frac{m_i}{q}}

where :math:`I_b` is the ion beam current.

The ideal specific impulse of an electrostatic thruster is:

.. math::

  I_{sp} &= \frac{v_e}{g_0}

         &= \frac{1}{g_0} \sqrt{2 V_b \frac{q}{m_i}}


As an example, use ``proptools`` to compute the thrust and specific impulse of singly charged Xenon ions with a beam voltage of 1000 V and a beam current of 1 A:

.. literalinclude:: examples/electric/thrust_isp.py

.. literalinclude:: examples/electric/thrust_isp_output.txt

This example illustrates the typical performance of electric propulsion systems: low thrust,  high specific impulse, and low mass flow.


Advantages over Chemical Propulsion
===================================
Electric propulsion is appealing because it enables higher specific impulse than chemical propulsion. The specific impulse of a chemical rocket depends on velocity to which a nozzle flow can be accelerated. This velocity is limited by the finite energy content of the working gas. In contrast, the particles leaving an electric thruster can be accelerated to very high velocities if sufficient electrical power is available.

In a chemical rocket, the kinetic energy of the exhaust gas is supplied by thermal energy released in a combustion reaction. For example, in the stoichiometric combustion of H\ :sub:`2` and O\ :sub:`2` releases an energy per particle of :math:`4.01 \times 10^{-19}` J, or 2.5 eV. If all of the released energy were converted into kinetic energy of the exhaust jet, the maximum possible exhaust velocity would be:

.. math::

  v_e = \sqrt{2 \frac{E}{m_{H_2O}}} = 5175 \, \mathrm{m s^{-1}}

corresponding to a maximum specific impulse of about 530 s. In practice, the most efficient flight engines have specific impulses of 450 to 500 s.

With electric propulsion, much higher energies per particle are possible. If we accelerate (singly) charged particles through a potential of 1000 V, the jet kinetic energy will be 1000 eV per ion. This is 400 times more energy per particle than is possible with chemical propulsion (or, for similar particle masses, 20 times higher specific impulse). Specific impulse in excess of 10,000 s is feasible with electric propulsion.


Power and Efficiency
====================

In an electric thruster, the energy to accelerate the particles in the jet must be supplied by an external source. Typically solar panels supply this electrical power, but some future concepts might use nuclear reactors. The available power limits the thrust and specific impulse of an electric thruster.

The kinetic power of the jet is given by:

.. math::

  P_{jet} = \frac{F^2}{2 \dot{m}}

where :math:`\dot{m}` is the jet mass flow. Increasing thrust with increasing mass flow (i.e. increasing :math:`I_{sp}`) will increase the kinetic power of the jet.

The power input required by the thruster (:math:`P_{in}`) is somewhat higher than the jet power. The ratio of the jet and input power is the total efficiency of the thruster, :math:`\eta_T`:

.. math::

  \eta_T \equiv \frac{P_{jet}}{P_{in}}

The total efficiency depends on several factors:

#. Thrust losses due to beam divergence. This loss is proportional to the cosine of the average beam divergence half-angle.
#. Thrust losses due to doubly charged ions in the beam. This loss is a function of the doubles-to-singles current ratio, :math:`I^{++}/I^+`
#. Thrust losses due to propellant gas escaping the thruster without being ionized. The fraction of the propellant mass flow which is ionized and enters the beam is the mass utilization efficiency, :math:`\eta_m`.
#. Electrical losses incurred in ion generation, power conversion, and powering auxiliary thruster components. These losses are captured by the electrical efficiency, :math:`\eta_e = \frac{I_b V_b}{P_{in}}`

Use ``proptools`` to compute the efficiency and required power of the example thruster. Assume that the beam divergence is :math:`\cos(10 \unicode{xb0})`, the double ion current fraction is 10%, the mass utilization efficiency is 90%, and the electrical efficiency is 85%:

.. literalinclude:: examples/electric/power.py

.. literalinclude:: examples/electric/power_output.txt

The overall efficiency of the thruster is about 70%. The required input power could be supplied by a few square meters of solar panels (at 1 AU from the sun).

The power, thrust, and specific impulse of a thruster are related by:

.. math::

  \frac{F}{P_{in}} = \frac{2 \eta_T}{g_0 I_{sp}}

Thus, for a power-constrained system the propulsion designer faces a trade-off between thrust and specific impulse.

.. plot:: examples/electric/plots/thrust_power_isp.py
  :include-source:
  :align: center


Optimal Specific Impulse
========================

For electrically propelled spacecraft, there is an optimal specific impulse which will maximize the payload mass fraction of a given mission. While increasing specific impulse decreases the required propellant mass, it also increases the required power at a particular thrust level, which increases the mass of the power supply. The optimal specific impulse minimizes the combined mass of the propellant and power supply.

The optimal specific impulse depends on several factors:

#. The mission thrust duration, :math:`t_m`. Longer thrust durations reduce the required thrust (if the same total impulse or :math:`\Delta v` is to be delivered), and therefore reduce the power and power supply mass at a given :math:`I_{sp}`. Therefore, longer thrust durations increase the optimal :math:`I_{sp}`.
#. The specific mass of the power supply, :math:`\alpha`. This is the ratio of power supply mass to power, and is typically 20 to 200 kg kW :sup:`-1` for solar-electric systems. The specific impulse optimization assumes that power supply mass is linear with respect to power. Increasing the specific mass reduces the optimal :math:`I_{sp}`.
#. The total efficiency of the thruster.
#. The :math:`\Delta v` of the mission. Higher :math:`\Delta v` (in a fixed time window) requires more thrust, and therefore leads to a lower optimal :math:`I_{sp}`.


Consider an example mission to circularize the orbit of a geostationary satellite launched onto a elliptical transfer orbit. Assume that the low-thrust circularization maneuver requires a :math:`\Delta v` of 2 km s :sup:`-1` over 100 days. The thruster is 70% efficient and the power supply specific mass is 50 kg kW :sup:`-1`:

.. literalinclude:: examples/electric/isp_opt.py

.. literalinclude:: examples/electric/isp_opt_output.txt


For the mathematical details of specific impulse optimization, see [Lozano]_.

.. [Lozano] P. Lozano, *16.522 Lecture Notes*, Lecture 3-4 Mission Analysis for Electric Propulsion. Online: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-522-space-propulsion-spring-2015/lecture-notes/MIT16_522S15_Lecture3-4.pdf
