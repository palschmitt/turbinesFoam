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
## @deftypefn {} {@var{retval} =} loading (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function loading(nfree,loads,order2)
    global nnodes ndof;
    global Ff Fr;
    totdof=nnodes*ndof;
    loadvector=zeros(1,totdof);
    icount=0;
    for inode=1:nnodes;
        for idof=1:ndof
            icount=icount+1;
            ivalue=loads(inode,idof);
            loadvector(icount)=ivalue;
        end
    end
    for idof=1:totdof
        ipos=order2(idof);
        if (ipos <= nfree)
            Ff(ipos)=loadvector(idof);
        else 
            Fr(ipos-nfree)=loadvector(idof);
        end
    end    
endfunction

