/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com

 * The most recent release of Palabos can be downloaded at
>>>>>>> bebbcca4e6a9e65b02d9b9bf9f4b48483d871410
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <mysrc/headersJo/buoyantProcessor3D.hh>
#include <mysrc/headersJo/buoyantProcessor3D_v2.hh>
#include <mysrc/headersJo/CouplingProcessor3D.hh>
#include <mysrc/headersJo/VolFracDensityCoupling.hh>
#include <mysrc/headersJo/TagProcessor3D.hh>
#include <mysrc/headersJo/ForceCoupling3D_FD.hh>
#include <mysrc/headersJo/DensVsedCoupling3D_FD.hh>
#include <mysrc/headersJo/Velocity_Coupling_FD.hh>
