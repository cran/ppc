C Output from Public domain Ratfor, version 1.0
      subroutine hclust1d(idebug,x,n,merge,height,order, val, val2)
      implicit double precision (a-h, o-z)
      integer idebug
      integer n,o
      integer merge(n-1,2),order(n)
      double precision x(n),height(n-1), val(n,3), val2(n,3)
      do23000 i=1,n
      order(i)=i
23000 continue
23001 continue
      do23002 i=1,n-1
      val(i,1)= -i
      val(i,2)=x(i+1)-x(i)
      val(i,3)=0
23002 continue
23003 continue
      val(n,1)= -n
      val(n,2)=99999
      val(n,3)=0
      lenval=n
      do23004 ii =1,(n - 1) 
      call whichism(val(1,2),lenval,o)
      height(ii) = val(o, 2)
      merge(ii,1)=val(o, 1)
      merge(ii,2)=val(o+1, 1)
      val(o, 3) = val(o, 2)
      val(o - 1, 2) = val(o - 1, 2) + val(o, 3)
      val(o, 2) = val(o + 1, 2) + val(o, 3)
      val(o, 1) = ii
      call mycopy1(val,val2,n,lenval,o+1)
      lenval=lenval-1
      call mycopy(val2,val,n,lenval)
      if(idebug.gt.0)then
      write(6,*) "merging",merge(ii,1),merge(ii,2)
      endif
23004 continue
23005 continue
      return
      end
      subroutine whichism(x, n,ans)
      implicit double precision (a-h, o-z)
      integer n, ans
      double precision x(n), minval
      minval= 10e9
      do23008 i = 1,n
      if(x(i).lt.minval)then
      ans=i
      minval=x(i)
      endif
23008 continue
23009 continue
      return
      end
      subroutine mycopy (x,x2,nrowx, n)
      implicit double precision (a-h, o-z)
      integer n,nrowx
      double precision x(nrowx,3) ,x2(nrowx,3)
      do23012 i=1,n
      do23014 j=1,3
      x2(i,j)=x(i,j)
23014 continue
23015 continue
23012 continue
23013 continue
      return
      end
      subroutine mycopy1 (x,x2,nrowx,n,ii)
      implicit double precision (a-h, o-z)
      integer n, ii, nrowx
      double precision x(nrowx,3) ,x2(nrowx,3)
      j=0
      do23016 i=1,n
      if(i.ne.ii)then
      j=j+1
      do23020 k=1,3
      x2(j,k)=x(i,k)
23020 continue
23021 continue
      endif
23016 continue
23017 continue
      return
      end
