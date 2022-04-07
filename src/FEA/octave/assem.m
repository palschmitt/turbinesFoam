## Copyright (C) 2020 Pal
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} assem (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function assem(Ke,elnodes)
    global K ndof nnode; 
    
    for innode=1:nnode
        ipos=(elnodes(innode)-1)*ndof;
        ipose=(innode-1)*ndof;
        for jnnode=1:nnode
            jpos=(elnodes(jnnode)-1)*ndof;
            jpose=(jnnode-1)*ndof;
            for idof=1:ndof
                for jdof=1:ndof
                    K(ipos+idof,jpos+jdof)=K(ipos+idof,jpos+jdof)+ ...
                        Ke(ipose+idof,jpose+jdof);
                end
            end
        end
    end
endfunction
