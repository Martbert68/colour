#include <time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BSIZE 28800000

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;
XImage *x_image;
unsigned char *x_buffer;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


void disp (unsigned char *,int,int);


/*void alloc_try (struct try *tri,int count)
{
	tri=(struct try *)malloc(sizeof (struct try)*count);
	int loop;
	for (loop=0;loop<count;loop++)
	{
		tri->pi=(float *)malloc(sizeof(float)*3);
		tri->pj=(float *)malloc(sizeof(float)*3);
		tri->pk=(float *)malloc(sizeof(float)*3);
	}
}*/

void fern(unsigned char *frame,int *origin,int iter,float len,float wig,float ang,float phi,int r, int g, int b)
{
	int op[2];
	float nlen;
	if (iter>1 && len >2)
	{
		iter --;	
		op[0]=origin[0]; op[1]=origin[1];
		linet( frame, op, (r*(600-len)/600) ,(g*(600-len)/600) , (b*(600-len)/600) , len , ang , 1+(len*len/2000) ,wig, len/10, 0);
		fern(frame,op,iter,len*0.7,wig,ang+0.4+phi,phi,r,g,b);
		op[0]=origin[0]; op[1]=origin[1];
		linet( frame, op, (r*(600-len)/600) ,(g*(600-len)/600) , (b*(600-len)/600) , len , ang , 1+(len*len/2000) ,wig, len/10, M_PI);
		fern(frame,op,iter,len*0.7,wig,ang-0.4+phi,phi,r,g,b);
	}
}



void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}

void fft (short *wav, float *mag, int freq, long length, long points)
{
	int p;
	for (p=0;p<points;p++)
	{
		long sel;
		float phi,st,ct;
		st=0;
		ct=0;
		for (sel=(length*p/points);sel<((length*(p+1))/points)-2;sel+=2)
		{
			phi=(float)freq*(float)sel*2*M_PI/(48000*2);	
			st+=(sin(phi)*(float)wav[sel]);
			ct+=(cos(phi)*(float)wav[sel]);
		}
		mag[p]=sqrt((st*st)+(ct*ct));
		if (p>0){ mag[p]+=mag[p-1];mag[p]/=2;}
		//printf ("P %d St %f \n",p,mag[p]);
	}
}

void get_amp(short *wav, float *amp,long length, long points)
{
	int p;
	for (p=0;p<points;p++)
	{
		long sel;
		long tot,count;
		tot=0;
		count=0;
		for (sel=(length*p/points);sel<((length*(p+1))/points)-2;sel+=2)
		{
			int val;
			val=wav[sel];
			count++;
			if (val>0){ tot+=val;}else{tot-=val;}
		}
		amp[p]=(float)tot/(9640*count);
		if (p>0){ amp[p]+=amp[p-1];amp[p]/=2;}
	}
}

int main(int argc,char *argv[])
{
	unsigned char *image1,*image2,*image3,*image4,*image5,*image6;
        image1=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image2=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image4=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image5=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // buffer
        image6=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // buffer

        long point,along;
	double m,n;
	init_x();

	int mins;
	int secs;
	int sec30;
	int rate;
	long bsize;
	int frate;
	int fcount;

	rate=48000;
	mins=10;
	frate=30;
	sec30=rate*2/30;
	

	secs=mins*60;
	bsize=secs*rate*2;
	fcount=0;

       	short *waveform;
        double *dwaveform;
	double chop;

        waveform=(short *)malloc(sizeof(short)*bsize);
        dwaveform=(double *)malloc(sizeof(double)*bsize);

	//load_image(image1,"scan0001.jpg",0);
	//load_image(image2,"scan0002.jpg",0);
	//load_image(image3,"scan0003.jpg",0);

	double llang,rrang,lang,rang,lphase,rphase,lf,rf,ld,rd,lfm,rfm,lfd,rfd,lfmf,rfmf;

	lf=157;
	rf=2*lf;

	lfd=5;
	rfd=5;

	ld=4;
	rd=7;

	lphase=-200*M_PI;rphase=200*M_PI;
	lang=0;rang=0;llang=0;rrang=0;
	clear (image4);


	//load_image(image2,"./back.jpg",1);
	//
	double  pl,pr,dl,dr,fil,fang,tt,mm,exp,bexp,bb,bbl,bbbl,bbexp,fl,ffl;

	fang=0;

	int trigger,mag,btrig,bbtrig;

	exp=1;
	bexp=1;
	bbexp=1;

	fl=100;
	ffl=100;

	bbtrig=0;


	clear (image4);

	for (point=0;point<bsize;point+=2)
	{
		m=(double)point/(double)bsize;
		n=(1-(2*m));
		if (n<0){n=-n;}
		//n=n*n*n;
		//n=n/2;
		double left,right,lv,rv,llv,rrv,ff;
		left=0;right=0;

		fang+=2*M_PI*((m))/48000;


		ff=540+(455*sin(n*m*2*M_PI*103));
		fil=10+(200*(1+(cos(fang))));

		lfmf=1911*n;
		rfmf=1810*n;

		lv=sin((n*lphase)+lang);
		rv=cos((n*rphase)+rang);

		if (lv>(1-n)){ lv=1;} if (lv<-(1-n)){ lv=-m;}
		if (rv>(1-n)){ rv=1;} if (rv<-(1-n)){ rv=-m;}

		llv=sin(llang);
		rrv=cos(rrang);

		if (llv>(1-n)){ llv=(m);} if (llv<-(1-n)){ llv=-(1);}
		if (rrv>(1-n)){ rrv=(m);} if (rrv<-(1-n)){ rrv=-(1);}

		lv+=(ff)*llv/1200;
		rv+=(ff)*rrv/1200;

		chop=cos(2*M_PI*365.25*m)*cos(m*2*M_PI*28)*cos(m*2*M_PI*7);
		chop=cos(2*M_PI*m*100*n);
		//chop=cos(2*M_PI*m*)*cos(m*2*M_PI*28);
		//chop=cos(m*2*M_PI*10);
		if (chop<0.5 && chop>0){chop=0;rf=ff*1.5;fang=0;}else{ffl=fl;fl=rf/3;}
		if (chop<0){chop=0;lf=ff;rf=ff*1.5;fang=0;}
		//lf=ff;rf=ff*1.5;
		//
		//

		left=10000*lv*chop;
		right=10000*rv*chop;


/*		dl=left-pl; dr=right-pr;

		if (dl>fil){left=pl+fil;} if (dl<-fil){left=pl-fil;}
		if (dr>fil){right=pr+fil;} if (dr<-fil){right=pr-fil;}


		pl=left;
		pr=right; */



		lang+=2*M_PI*(lf+(n*lfd*lfm)+(n*ld))/48000;
		rang+=2*M_PI*(rf-(n*rfd*rfm)-(n*rd))/48000;

		llang+=2*M_PI*(2*(lf+(n*n*lfd*rfm)))/48000;
		rrang+=2*M_PI*(2*(rf+(n*n*rfd*lfm)))/48000;

		lfm=sin(m*2*M_PI*lfmf);
		rfm=sin(m*2*M_PI*rfmf);

		double v,bl,br;
		int p1;


		if (point>96000)
		{
			v=0.8;
			p1=point-96000+2*((int)((m)*45000));

			bl=dwaveform[p1];
			br=dwaveform[p1+1];

			dl=bl-pl; dr=br-pr;

			if (dl>fil){bl=pl+fil;} if (dl<-fil){bl=pl-fil;}
			if (dr>fil){br=pr+fil;} if (dr<-fil){br=pr-fil;}

			pl=bl; pr=br;

			dwaveform[point]=left+(bl*v);
			dwaveform[point+1]=right+(br*v);


		}else{
			dwaveform[point]=left;
			dwaveform[point+1]=right;}



		waveform[point]=dwaveform[point];
		waveform[point+1]=dwaveform[point+1];


		//hat
		tt=cos(2*M_PI*m*600);
		bbbl=bbl;bbl=bb; bb=cos((M_PI/2)+2*M_PI*m*301)*(cos((3*M_PI/2)+2*M_PI*m*299));
		//if (bb>0.99 && bbtrig==0){bbtrig=1;bbexp=1;}
		//if (bb<0.99 && bbtrig==-1){bbtrig=0;}
		if (bb<bbl && bbbl<bbl) {bbtrig=1;bbexp=1;}
		if (tt>0.99 && trigger==0){trigger=1;exp=1;}
		if (tt<0.99 && trigger==-1){trigger=0;}
		if (tt<-0.99 && btrig==0){btrig=1;bexp=1;}
		if (tt>-0.99 && btrig==-1){btrig=0;}

		//waveform[point]=0;waveform[point+1]=0;chop=0;

		if (trigger>0 ){ 
			waveform[point]+=(rand()%12768*exp*(1)); waveform[point+1]+=(rand()%12768*exp*(1));
			trigger=2;
			if (exp<0.01){trigger=-1;}
			exp*=(0.999+(m*0.0009));
		}

		if (btrig>0 ){ 
			waveform[point]+=(rand()%12768*bexp*(1)); waveform[point+1]+=(rand()%12768*bexp*(1));
			btrig=2;
			if (bexp<0.01){btrig=-1;}
			bexp*=(0.9999-(m*0.0009));
		}

		//kick

		if (bbtrig>0 ){ 
			waveform[point]+=(1-chop)*22768*bbexp*sin(2*M_PI*(double)point*(60+(n*bb/30))/48000);
			waveform[point+1]+=(1-chop)*22768*bbexp*sin(2*M_PI*(double)point*(62-(m*bb/30))/48000); 
			if (bbexp<0.01){bbtrig=0;}
			bbexp*=(0.9995+(n*0.0003));
			bbtrig=2;
		}

		//drone
		waveform[point]+=bb*tt*(1-chop)*6768*sin(2*M_PI*(double)point*(ffl)/48000);
		waveform[point+1]+=bb*tt*(1-chop)*6768*sin(2*M_PI*(double)point*(ffl+4)/48000);

		if (point%sec30==0 && point>=sec30)
		//if (point%sec30==0 && point<-100)
		{
			//clear (image4);
			int ah,l,r,xx,xr,xl,lmv,rmv,yy,x,y;
			l=0;r=0;

			for (ah=0;ah<sec30;ah+=2)
			{
				lmv=waveform[point-sec30+ah];
				rmv=waveform[point-sec30+ah+1];
				if (lmv<0 && l==0){l=-1;}
				if (rmv<0 && r==0){r=-1;}
				if (l==-1 && lmv>0){l=ah;}
				if (r==-1 && rmv>0){r=ah;}
			}

			/*for (x=0;x<X_SIZE;x++)
			{
				for(y=0;y<Y_SIZE;y++)
				{
					plott(image4,x,y,rf*x/100,chop*y,tt*(x+y));
				}
			}*/

//void linet( unsigned char *image2, int *origin, int r, int g, int b, float length, float phi, float thick, float wig, float depth, float pha)
			int oo[2];
			oo[0]=(X_SIZE/2)+((n*X_SIZE/2)*sin(2*M_PI*m*20));
			oo[1]=(Y_SIZE/2)+((n*Y_SIZE/2)*cos(2*M_PI*m*20));
			linet(image4,oo,255,0,0,1500*bexp,tt*M_PI,8,lfm*10,30*chop,0);
			linet(image4,oo,255,0,0,1500*bexp,n*tt*M_PI/2,8,lfm*10,30*chop,0);
			oo[0]=(X_SIZE/2)+((n*X_SIZE/2)*sin(2*M_PI*m*18));
			oo[1]=(Y_SIZE/2)+((n*Y_SIZE/2)*cos(2*M_PI*m*18));
			linet(image4,oo,0,255,0,1500*bbexp,bb*M_PI,8,rfm*10,30*chop,0);
			linet(image4,oo,0,255,0,1500*bbexp,n*bb*M_PI/2,8,rfm*10,30*chop,0);
			oo[0]=(X_SIZE/2)+((n*X_SIZE/2)*sin(2*M_PI*m*30));
			oo[1]=(Y_SIZE/2)+((n*Y_SIZE/2)*cos(2*M_PI*m*15));
			linet(image4,oo,0,0,255,1500*exp,ff/200,8,fl/10,30*chop,0);
			linet(image4,oo,0,0,255,1500*exp,m*ff/200,8,fl/10,30*chop,0);


			merge(image4,image6,55);
			//ripple (unsigned char *painted, float xphase, float yphase, float xdepth, float ydepth, float xfreq, float yfreq)
			memcpy(image6,image4,X_SIZE*Y_SIZE*3);
			//merge(image4,image2,255*(n+10)/11);
			//ripple(image4,0,0,2,0,100*n,100*n);
			disp(image4,fcount++,1);
			if (chop){blurt(image4,20);}
		}
	}

        int *fhead,chan,sample_rate,bits_pers,byte_rate,ba,size;
        fhead=(int *)malloc(sizeof(int)*11);
        char fname[200];
        FILE *record;
        chan=2;
        sample_rate=48000;
        bits_pers=16;
        byte_rate=(sample_rate*chan*bits_pers)/8;
        ba=((chan*bits_pers)/8)+bits_pers*65536;

        fhead[0]=0x46464952;
        fhead[1]=36;
        fhead[2]=0x45564157;
        fhead[3]=0x20746d66;
        fhead[4]=16;
        fhead[5]=65536*chan+1;
        fhead[6]=sample_rate;
        fhead[7]=byte_rate;
        fhead[8]=ba;
        fhead[9]=0x61746164;
        fhead[10]=(bsize*chan*bits_pers)/8;


        record=fopen("crystal.wav","wb");
        fwrite(fhead,sizeof(int),11,record);
        fwrite(waveform,sizeof(short),bsize,record);
        fclose (record);

	close_x();

	exit(0);
}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*STRIDE;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/image%04d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

