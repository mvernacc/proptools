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

This example illustrates the typical performance of electric propulsion systems: low thrust and high specific impulse.


Advantages over Chemical Propulsion
==================================
Electric propulsion is appealing because it enables higher specific impulse than chemical propulsion. The specific impulse of a chemical rocket depends on velocity to which a nozzle flow can be accelerated. This velocity is limited by the finite energy content of the working gas. In contrast, the particles leaving an electric thruster can be accelerated to very high velocities if sufficient electrical power is available.

In a chemical rocket, the kinetic energy of the exhaust gas is supplied by thermal energy released in a combustion reaction. For example, in the stoichiometric combustion of H\ :sub:`2` and O\ :sub:`2` releases an energy per particle of :math:`4.01 \times 10^{-19}` J, or 2.5 eV. If all of the released energy were converted into kinetic energy of the exhaust jet, the maximum possible exhaust velocity would be:

.. math::

  v_e = \sqrt{2 \frac{E}{m_{H_2O}}} = 5175 \, \mathrm{m s^{-1}}

corresponding to a maximum specific impulse of about 530 s. In practice, the most efficient flight engines have specific impulses of 450 to 500 s.

With electric propulsion, much higher energies per particle are possible. If we accelerate (singly) charged particles through a potential of 1000 V, the jet kinetic energy will be 1000 eV per ion. This is 400 times more energy per particle than is possible with chemical propulsion (or, for similar particle masses, 20 times higher specific impulse). Specific impulse in excess of 10,000 s is feasible with electric propulsion.


Power
=====



Optimal Specific Impulse
========================
