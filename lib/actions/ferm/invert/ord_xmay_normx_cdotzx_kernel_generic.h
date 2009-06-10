void ord_xmay_normx_cdotzx_kernel(int lo, int hi, int my_id, ord_xmay_normx_cdotzx_arg *a) 
{
  
  REAL32* x_ptr=&(a->x_ptr[lo]);
  REAL32* y_ptr=&(a->y_ptr[lo]);
  REAL32* z_ptr=&(a->z_ptr[lo]);
  REAL32 a_re = a->a_re;
  REAL32 a_im = a->a_im;
  REAL64 norm_array[3]={0,0,0};
  
  int len = hi-lo;
  for(int count = 0; count < len; count+=4) { 

    x_ptr[count] -= a_re*y_ptr[count];
    x_ptr[count] += a_im*y_ptr[count+1];
    
    x_ptr[count+1] -= a_im*y_ptr[count];
    x_ptr[count+1] -= a_re*y_ptr[count+1];
    
    x_ptr[count+2] -= a_re*y_ptr[count+2];
    x_ptr[count+2] += a_im*y_ptr[count+3];
    
    x_ptr[count+3] -= a_im*y_ptr[count+2];
    x_ptr[count+3] -= a_re*y_ptr[count+3];
    
    norm_array[0] += x_ptr[count]*x_ptr[count];
    norm_array[0] += x_ptr[count+1]*x_ptr[count+1];
    norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
    norm_array[0] += x_ptr[count+3]*x_ptr[count+3];
    
    norm_array[1] += z_ptr[count]*x_ptr[count];
    norm_array[1] += z_ptr[count+1]*x_ptr[count+1];
    norm_array[1] += z_ptr[count+2]*x_ptr[count+2];
    norm_array[1] += z_ptr[count+3]*x_ptr[count+3];
    
    norm_array[2] += z_ptr[count]*x_ptr[count+1];
    norm_array[2] -= z_ptr[count+1]*x_ptr[count];
    norm_array[2] += z_ptr[count+2]*x_ptr[count+3];
    norm_array[2] -= z_ptr[count+3]*x_ptr[count+2];
    
  }
  a->norm_space[3*my_id]=norm_array[0];
  a->norm_space[3*my_id+1]=norm_array[1];
  a->norm_space[3*my_id+2]=norm_array[2];
}      
