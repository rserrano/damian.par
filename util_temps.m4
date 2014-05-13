define(`dimc',`ifelse($1,`1',`:',`:,dimc(decr($1))')')
define(`dim',`dimc($1)')

define(`name_func',`$1$2$3')
define(`name_concat',`name_func(`$1',substr(`$2',`0',`1'),$3)')

define(`shlistc',`ifelse($3,$1,`$2($3)',`$2($3),shlistc($1,$2,incr($3))')')
define(`shlist',`shlistc($1,$2,1)')

define(`ctofort_func', `subroutine name_concat(`ctofort',$1,$2)(caddr, array, sh)
type(c_ptr), intent(in) :: caddr
integer, dimension(:), intent(in) :: sh
$1, dimension(dim($2)), allocatable, intent(out) :: array
$1, pointer, dimension(dim($2)) :: fptr
call c_f_pointer(caddr, fptr, sh)
allocate(array(shlist($2,sh)))
array = fptr
call c_free(caddr)
end subroutine
')

ctofort_func(`real*8', `1')
ctofort_func(`real*8', `2')
ctofort_func(`integer', `1')
ctofort_func(`integer', `2')
ctofort_func(`character*1', `1')
ctofort_func(`character*1', `2')
ctofort_func(`logical*1', `1')
ctofort_func(`logical*1', `2')

