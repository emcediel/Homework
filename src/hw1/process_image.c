#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    x = x>=im.w ? (im.w-1) : x;
    x = x<0 ? 0 : x;
    y = y>=im.h ? (im.h-1) : y;
    y = y<0 ? 0 : y;
    return im.data[x+(y*im.w)+(c*im.w*im.h)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if (x<im.w && y<im.h && x>=0 & y>=0){
        im.data[x+(y*im.w)+(c*im.w*im.h)]=v;
    }
    
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data,im.data,im.w*im.h*im.c*sizeof(int));
    // TODO Fill this in
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
            float R=get_pixel(im,i,j,0);
            float G=get_pixel(im,i,j,1);
            float B=get_pixel(im,i,j,2);
            gray.data[i+j*im.w]=0.299*R + 0.587*G + .114*B;
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
            set_pixel(im,i,j,c,get_pixel(im,i,j,c)+v);
        }
    }
    
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i=0;i<im.w;i++)
    {
        for (int j=0;j<im.h;j++)
        {
            for (int c=0;c<im.c;c++)
            {
                if (get_pixel(im,i,j,c)<0)
                {
                    set_pixel(im,i,j,c,0);
                } 
                else if (get_pixel(im,i,j,c)>1)
                {
                    set_pixel(im,i,j,c,1);
                }
            }   
        }
        
    }
    
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
            float R = get_pixel(im,i,j,0);
            float G = get_pixel(im,i,j,1);
            float B = get_pixel(im,i,j,2);
            float V = three_way_max(R,G,B); //We calculate Value
            float C = V - three_way_min(R,G,B); 
            float S = V==0 ? 0 : C/V; //We calculate Saturation
            float H_prime = C==0 ? 0 : (V==R ? (G-B)/C :(V==G ? ((B-R)/C) +2 : ((R-G)/C)+4)); 
            float H = H_prime < 0 ? (H_prime/6)+1 : (H_prime/6);   //We calculate Hue
            set_pixel(im,i,j,0,H);
            set_pixel(im,i,j,1,S);
            set_pixel(im,i,j,2,V);
        }
    }
    
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int i=0;i<=im.w-1;i++)
    {
        for(int j=0;j<=im.h-1;j++)
        {
            float H = get_pixel(im,i,j,0)*6;
            float V = get_pixel(im,i,j,2);
            float S = get_pixel(im,i,j,1);
            float Hi = floor(H);
            float F = H - Hi;
            float P = V*(1-S);
            float Q = V*(1-F*S);
            float T = V*(1-(1-F)*S);
            float R;
            float G;
            float B;
            if (Hi == 0)
            {
                R=V;
                G=T;
                B=P;
            }
            else if(Hi==1)
            {
                R=Q;
                G=V;
                B=P;
            }
            else if(Hi==2)
            {
                R=P;
                G=V;
                B=T;
            }
            else if(Hi==3)
            {
                R=P;
                G=Q;
                B=V;
            }
            else if(Hi==4)
            {
                R=T;
                G=P;
                B=V;
            }
            else
            {
                R=V;
                G=P;
                B=Q;
            }
            set_pixel(im,i,j,0,R);
            set_pixel(im,i,j,1,G);
            set_pixel(im,i,j,2,B);

        }
    }

}
void scale_image(image im, int c, float v)
{
    for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
            set_pixel(im,i,j,c,get_pixel(im,i,j,c)*v);
        }
    }
}