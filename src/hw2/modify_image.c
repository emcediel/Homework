#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    int nn_x=round(x);
    int nn_y=round(y);
    float v=get_pixel(im,nn_x,nn_y,c);
    return v;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image new_im = make_image(w, h, im.c);
    float aX= (float)im.w/ (float)w;
    float aY= (float)im.h/ (float)h;
    float bX=0.5*aX-0.5;
    float bY=0.5*aY-0.5; 
    for (int i=0;i<w;i++)
    {
        for(int j=0;j<h;j++)
        {
          for(int k=0;k<im.c;k++)
          {
            float srcX=aX*i+bX;
            float srcY=aY*j+bY;
            float interVal=nn_interpolate(im,srcX,srcY,k);
            set_pixel(new_im,i,j,k,interVal);
          }
        }
    }
    return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    int f_x=floor(x);
    int f_y=floor(y);
    int c_x=ceil(x);
    int c_y=ceil(y);
    float q1=get_pixel(im,f_x,f_y,c);
    float q2=get_pixel(im,c_x,f_y,c);
    float q3=get_pixel(im,c_x,c_y,c);
    float q4=get_pixel(im,f_x,c_y,c);
    float r1=(x-f_x)*(y-f_y);
    float r2=(c_x-x)*(y-f_y);
    float r3=(c_x-x)*(c_y-y);
    float r4=(x-f_x)*(c_y-y);
    float v=q1*r3+q2*r4+q3*r1+q4*r2;
    return v;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
   image new_im = make_image(w, h, im.c);
    float aX= (float)im.w/ (float)w;
    float aY= (float)im.h/ (float)h;
    float bX=0.5*aX-0.5;
    float bY=0.5*aY-0.5; 
    for (int i=0;i<w;i++)
    {
        for(int j=0;j<h;j++)
        {
          for(int k=0;k<im.c;k++)
          {
            float srcX=aX*i+bX;
            float srcY=aY*j+bY;
            float interVal=bilinear_interpolate(im,srcX,srcY,k);
            set_pixel(new_im,i,j,k,interVal);
          }
        }
    }
    return new_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
   float sum=0;
   for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
          for (int k=0;k<im.c;k++)
          {
              sum+=get_pixel(im,i,j,k);
          }
          
      }
      
    }
    for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
          for (int k=0;k<im.c;k++)
          {
              set_pixel(im,i,j,k,(float) get_pixel(im,i,j,k)/sum);
          }
          
      }
      
    }

}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image new_im = make_image(w, w, 1);
    
    for (int j=0;j<new_im.h;j++)
        {
          for (int i=0;i<new_im.w;i++)
          {
              set_pixel(new_im,i,j,0,1);
          }
          
        }
    
    
   l1_normalize(new_im);
    
    return new_im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
   assert(im.c==filter.c || filter.c==1);
   float val=0;
   float f_pixel, im_pixel;
   image new_im;
   int offsetW,offsetH;
   if (im.c==filter.c && preserve == 0)
   {
    new_im = make_image(im.w, im.h, 1);
     for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
          val=0;
          for(int k=0;k<im.c;k++)
          {
            for (int l=0;l<filter.h;l++)
            {
              for (int w=0;w<filter.w;w++)
              {
                  f_pixel = get_pixel(filter,w,l,k);
                  offsetW = (int) i-filter.w/2 + w;
                  offsetH = (int) j-filter.h/2 +l;
                  im_pixel = get_pixel(im,offsetW,offsetH,k);
                  val += f_pixel * im_pixel;
              }          
            }
          }
          set_pixel(new_im,i,j,0, val);
        }
    }
   }
   else if (im.c==filter.c && preserve == 1)
   {
     new_im = make_image(im.w, im.h, im.c);
     for (int k=0;k<im.c;k++)
      {
        for(int j=0;j<im.h;j++)
        {
          for(int i=0;i<im.w;i++)
          {
            val=0;
            for (int l=0;l<filter.h;l++)
            {
              for (int w=0;w<filter.w;w++)
              {
                  f_pixel = get_pixel(filter,w,l,k);
                  offsetW = (int) i-filter.w/2 + w;
                  offsetH = (int) j-filter.h/2 +l;
                  im_pixel = get_pixel(im,offsetW,offsetH,k);
                  val += f_pixel * im_pixel;
              }
            }
            set_pixel(new_im,i,j,k,(float)val);
          }
        }
        
      }
   }
  else if (im.c!=filter.c && preserve == 0)
   {//im.c != filter.c && preserve == 0,put eveything into one channel
       new_im = make_image(im.w, im.h, 1);
        for (int i=0;i<im.w;i++)
        {
            for(int j=0;j<im.h;j++)
            {
              val=0;
              for(int k=0;k<im.c;k++)
              {
                for (int l=0;l<filter.h;l++)
                {
                  for (int w=0;w<filter.w;w++)
                  {
                  f_pixel = get_pixel(filter,w,l,0);
                  offsetW = (int) i-filter.w/2 + w;
                  offsetH = (int) j-filter.h/2 +l;
                  im_pixel = get_pixel(im,offsetW,offsetH,k);
                  val += f_pixel * im_pixel;
                  }
                }
              }
              set_pixel(new_im,i,j,0,val);
            }
          }
   }
     else
    {
     new_im = make_image(im.w, im.h, im.c);
      for (int k=0;k<im.c;k++)
        {
          for(int j=0;j<im.h;j++)
          {
            for(int i=0;i<im.w;i++)
            {
              val=0;
              for (int l=0;l<filter.h; l++)
              {
                for (int w=0;w<filter.w; w++)
                {
                  f_pixel = get_pixel(filter,w,l,0);
                  offsetW = (int) i-filter.w/2 + w;
                  offsetH = (int) j-filter.h/2 +l;
                  im_pixel = get_pixel(im,offsetW,offsetH,k);
                  val += f_pixel * im_pixel;
                }       
              }
              set_pixel(new_im,i,j,k,val);
            }
          }
        }  
    }
    return new_im;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
   image new_im = make_image(3, 3, 1);
   int filter[9]={0,-1,0,-1,4,-1,0,-1,0};
   for (int i=0;i<(new_im.h*new_im.w);i++)
   {
     new_im.data[i] = filter[i];
   }
    return new_im;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
   image new_im = make_image(3, 3, 1);
   int filter[9]={0,-1,0,-1,5,-1,0,-1,0};
   for (int i=0;i<(new_im.h*new_im.w);i++)
   {
     new_im.data[i] = filter[i];
   }
    return new_im;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
   image new_im = make_image(3, 3, 1);
   int filter[9]={-2,-1,0,-1,1,1,0,1,2};
   for (int i=0;i<(new_im.h*new_im.w);i++)
   {
     new_im.data[i] = filter[i];
   }
    return new_im;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: 
// We should use preserve on Sharpen, Box and Emboss filters given that these are kernels that we use to enhance an image, 
//hence it is more for a user that would like to sharpen or style their images, so they would like to have a rgb version 
//of it so they don't loose color while applying these filters.  On the other hand, Highpass filter is mainly used to find
//edges, in this scenario we don't care about the colors of the image but only about the edges, and this might be used as 
//a preprocessing step for feature extraction or edge detection, so rgb will not be necessary or it would even be confusing 
//as we would have identified edges on each channel (3 answers) instead of only one general answer in 1 channel



// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer:
// We should do post-processing to filters that are changing the overall magnitude of our image, this means that filters that have 
//large numbers outside the [0,1] range will tend to change the overall magnitude of the image and post-processing such as normalization will be required
//hence, we would recommend post-processing for Sharpen, Emboss and Highpass filters (also for Sobel filters that we will do later). On the other hand
//post-processing will not be necessary for box filters as they are between the range [0,1]. 


image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
   int size= ceil(sigma*6) + ((int)ceil(sigma*6) % 2 ^ 1);
   
   int offsetW,offsetH;  
   float gauss;
   image new_im = make_image(size, size, 1);
   for (int j=0;j<new_im.h;j++)
   {
     for(int i=0;i<new_im.w;i++)
     {
        offsetW = (int) i-new_im.w/2;
        offsetH = (int) j-new_im.h/2;
        gauss = (1/(2*M_PI*sigma*sigma))*exp(-1*(((offsetW*offsetW)+(offsetH*offsetH))/(2*sigma*sigma)));
        set_pixel(new_im,i,j,0,gauss);
     }
     
   }


    l1_normalize(new_im);
    return new_im;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
   assert(a.w == b.w && a.h == b.h && a.c == b.c);
   image new_im = make_image(a.w, a.h, a.c);
   for (int i=0;i<a.w;i++)
    {
      for (int j=0;j<a.h;j++)
      {
          for (int k=0;k<a.c;k++)
          {
              set_pixel(new_im,i,j,k,(float) (get_pixel(a,i,j,k)+get_pixel(b,i,j,k)));
          }
          
      }
      
    }
    return new_im;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
   assert(a.w == b.w && a.h == b.h && a.c == b.c);
   image new_im = make_image(a.w, a.h, a.c);
   for (int i=0;i<a.w;i++)
    {
      for (int j=0;j<a.h;j++)
      {
          for (int k=0;k<a.c;k++)
          {
              set_pixel(new_im,i,j,k,(float) (get_pixel(a,i,j,k)-get_pixel(b,i,j,k)));
          }
          
      }
      
    }
    return new_im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
   image new_im = make_image(3, 3, 1);
   int filter[9]={-1,0,1,-2,0,2,-1,0,1};
   for (int i=0;i<(new_im.h*new_im.w);i++)
   {
     new_im.data[i] = filter[i];
   }
    return new_im;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image new_im = make_image(3, 3, 1);
    int filter[9]={-1,-2,-1,0,0,0,1,2,1};
    for (int i=0;i<(new_im.h*new_im.w);i++)
    {
      new_im.data[i] = filter[i];
    }
      return new_im;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
   float min=get_pixel(im,0,0,0),max=get_pixel(im,0,0,0);

   for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
          for (int k=0;k<im.c;k++)
          {
              if (get_pixel(im,i,j,k)>max)
              {
                max=get_pixel(im,i,j,k);
              }
              if (get_pixel(im,i,j,k)<min)
              {
                min=get_pixel(im,i,j,k);
              }
          }
          
      }
      
    }

    for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
          for (int k=0;k<im.c;k++)
          {
              set_pixel(im,i,j,k,(float)(get_pixel(im,i,j,k)-min)/(max-min));
          }
          }
          
      }
      
    }




image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image mag = make_image(im.w, im.h, 1);
    image grad = make_image(im.w, im.h, 1);
    image s_x=make_gx_filter();
    image s_y=make_gy_filter();
    image G_x=convolve_image(im, s_x, 0);
    image G_y=convolve_image(im, s_y, 0);
    float magnitude,gradient;
    for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
        magnitude=sqrt(get_pixel(G_x,i,j,0)*get_pixel(G_x,i,j,0)+get_pixel(G_y,i,j,0)*get_pixel(G_y,i,j,0));
        gradient=atan2(get_pixel(G_y,i,j,0),get_pixel(G_x,i,j,0));
        set_pixel(mag,i,j,0,magnitude);
        set_pixel(grad,i,j,0,gradient);
      }
          
    }
    image *sobelimg = calloc(2, sizeof(image));
    sobelimg[0]=mag;
    sobelimg[1]=grad;
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  image f = make_gaussian_filter(4);
  im = convolve_image(im, f, 1);
  image *sobelimg = sobel_image(im);
  image mag=sobelimg[0];
  feature_normalize(mag);
  image grad=sobelimg[1];
  feature_normalize(grad);
  image new_im = make_image(im.w, im.h, 3);
  for (int i=0;i<im.w;i++)
    {
      for (int j=0;j<im.h;j++)
      {
        set_pixel(new_im,i,j,0,(float)(get_pixel(grad,i,j,0))); //we map gradient to hue (c=0)
        set_pixel(new_im,i,j,1,(float)(get_pixel(mag,i,j,0))); //we map magnitude to saturation (c=1)
        set_pixel(new_im,i,j,2,(float)(get_pixel(mag,i,j,0))); //we map value to value (c=2)
      }
          
      }
      
  hsv_to_rgb(new_im);

  return new_im;
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/
