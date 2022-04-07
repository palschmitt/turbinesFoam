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
## @deftypefn {} {@var{retval} =} neworder (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function [nfree,order1,order2]=neworder(restraints)
    global nnodes ndof;
    
    totdof=nnodes*ndof;
    order1=zeros(1,totdof);
    order2=zeros(1,totdof);
    irestraint=0;
    icount=0;
    restraintlist=zeros(1,totdof);
    for inode=1:nnodes;
        for idof=1:ndof
            icount=icount+1;
            ivalue=restraints(inode,idof);
            restraintlist(icount)=ivalue;
            if(ivalue == 1)
                irestraint=irestraint+1;
            end
        end
    end
    nrestraints=irestraint;
    nfree=totdof-nrestraints;
   
    irestraint=0;
    ifree=0;
    for idof=1:totdof;
        if(restraintlist(idof) == 1)
            irestraint=irestraint+1;
            order1(nfree+irestraint)=idof;
            order2(idof)=nfree+irestraint;
        else
            ifree=ifree+1;
            order1(ifree)=idof;
            order2(idof)=ifree;
        end
    end   
endfunction
