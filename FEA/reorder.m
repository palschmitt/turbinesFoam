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
## @deftypefn {} {@var{retval} =} reorder (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function reorder(nfree,order2)
global K Kff Kfr Krf Krr;
global nnodes ndof;
totdof=nnodes*ndof;
    for idof=1:totdof
        for jdof=1:totdof
            ipos=order2(idof);
            jpos=order2(jdof);
            if (ipos <= nfree && jpos <=nfree)
                Kff(ipos,jpos)=K(idof,jdof);
            elseif (ipos > nfree && jpos <=nfree)
                Krf(ipos-nfree,jpos)=K(idof,jdof);
            elseif (ipos <= nfree && jpos > nfree)
                Kfr(ipos,jpos-nfree)=K(idof,jdof);
            else
                Krr(ipos-nfree,jpos-nfree)=K(idof,jdof);
            end
        end
    end

endfunction
