/*
 * xmippImagic.cc - source file for ImageImagicT<T> class
 */

#include "../xmippImagic.hh"

/* Constants for the offset into the Imagic header file */
static const unsigned IMAGIC_IFOL_OFFSET=4, IMAGIC_IXLP_OFFSET=48,
  IMAGIC_IYLP_OFFSET=52, IMAGIC_TYPE_OFFSET=56, IMAGIC_WORD_LEN=4,
  IMAGIC_RECORD_LEN=1024;

/* Constants defining the Imagic header and some of its fields */
static const unsigned int IMAGIC_HEADER_BLOCK_SIZE = 256;
static const unsigned int IMAGIC_IDX_IMN = 0,
  IMAGIC_IDX_IFOL = 1, IMAGIC_IDX_NHFR = 3,
  IMAGIC_IDX_NDATE = 5, IMAGIC_IDX_NMONTH = 4, // These seem reversed in the spec
  IMAGIC_IDX_NYEAR = 6, IMAGIC_IDX_NHOUR = 7,
  IMAGIC_IDX_NMINUT = 8, IMAGIC_IDX_NSEC = 9,
  IMAGIC_IDX_NPIX2 = 10, IMAGIC_IDX_NPIXEL = 11,
  IMAGIC_IDX_IXLP1 = 12, IMAGIC_IDX_IYLP1 = 13,
  IMAGIC_IDX_TYPE = 14, IMAGIC_IDX_NAME = 29, IMAGIC_IDX_NAMELEN = 80,
  IMAGIC_IDX_ARCHTYPE = 68;


template <class T>
void ImageImagicT<T>::parseFname()
{
  if (!name_parsed)
  {
    hedfname = "";
    imgfname = "";
    imgnum = -1;
    // Look for special IMAGIC format: 'imagic:<hedfile>:<imgnum>'
    if (fn_img.find (IMAGIC_TAG) != string::npos)
    {
      const string::size_type imgnumpos = fn_img.rfind (IMAGIC_TAG_SEP);
      if (imgnumpos > IMAGIC_TAG_LEN)
      {
	hedfname = fn_img.substr (IMAGIC_TAG_LEN, imgnumpos-IMAGIC_TAG_LEN);
	imgfname = hedfname.substitute_extension (IMAGIC_HEADER_EXT,
						  IMAGIC_IMAGE_EXT);
	imgnum = atoi ((fn_img.substr (imgnumpos+1)).c_str());
      }
    }
    name_parsed = true;
  }
}

template <class T>
bool ImageImagicT<T>::read (const FileName &name) _THROW
{
  rename (name);
  ImageImagicInfo img_info = ImagicGetImgInfo (getHedFname());

  FileName img_fname = getImgFname();
  if (img_fname == "")
    REPORT_ERROR (1501, "ImageImagic::read: File " + name +
		  " doesn't seem fit Imagic format");
  FILE *img_fh;
  if ((img_fh = fopen (img_fname.c_str(), "rb")) == NULL)
    REPORT_ERROR (1501,"ImageImagic::read: IMAGIC file " + img_fname +
		  " not found");

  const int imgnum = getImgNum();
  const size_t img_offset = img_info.xsize * img_info.ysize;
  // Read the image data
  img.resize (img_info.ysize, img_info.xsize);
  const bool reversed = false;
  switch (img_info.img_types[imgnum])
  {
  case IMAGIC_REAL:
    {
      float data;
      const unsigned size = 4;
      fseek (img_fh, imgnum*size*img_offset, SEEK_SET);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
	{
	  FREAD (&data, size, 1, img_fh, reversed);
	  MULTIDIM_ELEM(img,i) = data;
	}
      break;
    }
  case IMAGIC_INTG:
    {
      short int data;
      const unsigned size = 2;
      fseek (img_fh, imgnum*size*img_offset, SEEK_SET);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
	{
	  FREAD (&data, size, 1, img_fh, reversed);
	  MULTIDIM_ELEM(img,i) = data;
	}
      break;
    }
  case IMAGIC_PACK:
    {
      unsigned char data;
      const unsigned size = 1;
      fseek (img_fh, imgnum*size*img_offset, SEEK_SET);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
	{
	  FREAD (&data, size, 1, img_fh, reversed);
	  MULTIDIM_ELEM(img,i) = data;
	}
      break;
    }
  default:
    REPORT_ERROR(1501,"ImageImagicType not supported for this imgtype!");
    break;
  }
  fclose (img_fh);
  return (TRUE);
}

bool ImageImagicT<complex<double> >::read (const FileName &name) _THROW
{
  rename (name);
  ImageImagicInfo img_info = ImagicGetImgInfo (getHedFname());

  FileName img_fname = getImgFname();
  if (img_fname == "")
    REPORT_ERROR (1501, "ImageImagic::read: File " + name +
		  " doesn't seem fit Imagic format");
  FILE *img_fh;
  if ((img_fh = fopen (img_fname.c_str(), "rb")) == NULL)
    REPORT_ERROR (1501,"ImageImagic::read: IMAGIC file " + img_fname +
		  " not found");

  const int imgnum = getImgNum();
  const size_t img_offset = img_info.xsize * img_info.ysize;
  // Read the image data
  const bool reversed = false;
  switch (img_info.img_types[imgnum])
  {
  case IMAGIC_COMP:
    {
      float a, b;
      const unsigned size = 4;
      fseek (img_fh, imgnum*(size*2)*img_offset, SEEK_SET);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
	{
	  // read real part of a complex number
	  FREAD (&a, size, 1, img_fh, reversed);
	  // read imaginary part of a complex number
	  FREAD (&b, size, 1, img_fh, reversed);
	  // Assign the number
	  complex<double> c(a,b);
	  MULTIDIM_ELEM(img,i)=c;
	}
      break;
    }
  default:
    REPORT_ERROR(1501,"ImageImagicType not supported for this imgtype!");
    break;
  }
  fclose (img_fh);
  return (TRUE);
}

template <class T>
bool ImagicWriteImagicFile (const FileName &hed_fname,
			    const vector<ImageT<T> *> &imgs,
			    ImageImagicType img_type)
{
  const FileName img_fname = hed_fname.substitute_extension (IMAGIC_HEADER_EXT,
							     IMAGIC_IMAGE_EXT);
  FILE *imagic_hed = fopen (hed_fname.c_str(), "wb");
  FILE *imagic_img = fopen (img_fname.c_str(), "wb");
  if (imagic_hed && imagic_img)
  {
    // Write the header information
    int header_block[IMAGIC_HEADER_BLOCK_SIZE];
    for (unsigned int imgcount = 0; imgcount < imgs.size(); imgcount++)
    {
      const ImageT<T> *image = imgs[imgcount];
      if (!image)
      {
	if (imagic_hed) fclose (imagic_hed);
	if (imagic_img) fclose (imagic_img);
	return (FALSE);
      }
      memset (header_block, 0, sizeof (header_block));
      header_block[IMAGIC_IDX_IMN] = imgcount+1;
      header_block[IMAGIC_IDX_IFOL] = imgs.size()-(imgcount+1);
      header_block[IMAGIC_IDX_NHFR] = 1;
      const time_t nowt = time(NULL);
      const struct tm *nowtm = localtime (&nowt);
      header_block[IMAGIC_IDX_NMONTH] = nowtm->tm_mon+1;
      header_block[IMAGIC_IDX_NDATE] = nowtm->tm_mday;
      header_block[IMAGIC_IDX_NYEAR] = nowtm->tm_year+1900;
      header_block[IMAGIC_IDX_NHOUR] = nowtm->tm_hour;
      header_block[IMAGIC_IDX_NMINUT] = nowtm->tm_min;
      header_block[IMAGIC_IDX_NSEC] = nowtm->tm_sec;
      header_block[IMAGIC_IDX_NPIX2] = XSIZE((*image)())*YSIZE((*image)());
      header_block[IMAGIC_IDX_NPIXEL] = header_block[IMAGIC_IDX_NPIX2];
      header_block[IMAGIC_IDX_IXLP1] = XSIZE((*image)());
      header_block[IMAGIC_IDX_IYLP1] = YSIZE((*image)());
      string formatstr;
      switch (img_type)
      {
      case IMAGIC_REAL:
	formatstr = "REAL";
	break;
      case IMAGIC_INTG:
	formatstr = "INTG";
	break;
      case IMAGIC_PACK:
	formatstr = "PACK";
	break;
      case IMAGIC_COMP:
	formatstr = "COMP";
	break;
      }
      memcpy (&header_block[IMAGIC_IDX_TYPE], formatstr.c_str(), 4);
      strncpy ((char *) &header_block[IMAGIC_IDX_NAME],
	       image->name().c_str(), IMAGIC_IDX_NAMELEN);
#if defined(_LINUX)
      static const unsigned int ARCH_VAL = 33686018;
#elif defined(_SUN)
      static const unsigned int ARCH_VAL = 67372036;
#else
      static const unsigned int ARCH_VAL = 33686018;
#endif
      // This next line will generate an error if not using linux or sun!
      header_block[IMAGIC_IDX_ARCHTYPE] = ARCH_VAL;

      fwrite (header_block, sizeof (header_block), 1, imagic_hed);
      // Write the image data to the .img file
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY((*image)())
	{
	  switch (img_type)
	  {
	  case IMAGIC_REAL:
	    {
	      const float p = (float) MULTIDIM_ELEM ((*image)(), i);
	      FWRITE (&p, sizeof (p), 1, imagic_img, FALSE);
	      break;
	    }
	  case IMAGIC_INTG:
	    {
	      const unsigned short p =
		(unsigned short) MULTIDIM_ELEM ((*image)(), i);
	      FWRITE (&p, sizeof (p), 1, imagic_img, FALSE);
	      break;
	    }
	  case IMAGIC_PACK:
	  case IMAGIC_COMP:
	    // NT: TODO: implement these
	    fprintf (stderr, "Unsupported format for writeSelFile!\n");
	    fclose (imagic_hed);
	    fclose (imagic_img);
	    return (FALSE);
	    break;
	  }
	}
    }
  }
  if (imagic_hed) fclose (imagic_hed);
  if (imagic_img) fclose (imagic_img);
  return (TRUE);
}


const ImageImagicInfo ImagicGetImgInfo (const FileName &hed_fname)
{
  ImageImagicInfo info;

  FILE *fp;
  if ((fp = fopen (hed_fname.c_str(), "rb")) != NULL)
  {
    // Read how many images (IFOL) and img size (NPIXEL)
    fseek (fp, IMAGIC_IFOL_OFFSET, SEEK_SET);
    unsigned int num_img;
    fread (&num_img, IMAGIC_WORD_LEN, 1, fp);
    num_img++; // Don't forget to include the current image
    info.num_img = num_img;
    fseek (fp, IMAGIC_IXLP_OFFSET, SEEK_SET);
    fread (&info.xsize, IMAGIC_WORD_LEN, 1, fp);
    fseek (fp, IMAGIC_IYLP_OFFSET, SEEK_SET);
    fread (&info.ysize, IMAGIC_WORD_LEN, 1, fp);
    for (unsigned int i = 0; i < num_img; i++)
    {
      // Determine what data type for this image
      char typeval[IMAGIC_WORD_LEN+1];
      fseek (fp, i*IMAGIC_RECORD_LEN+IMAGIC_TYPE_OFFSET, SEEK_SET);
      fread (&typeval, IMAGIC_WORD_LEN, 1, fp);
      typeval[IMAGIC_WORD_LEN]='\0';
      // fIform values are defined in headerXmipp::set_header
      if (strncmp (typeval, "REAL", 4) == 0)
	info.img_types.push_back (IMAGIC_REAL);
      else if (strncmp (typeval, "INTG", 4) == 0)
	info.img_types.push_back (IMAGIC_INTG);
      else if (strncmp (typeval, "PACK", 4) == 0)
	info.img_types.push_back (IMAGIC_PACK);
      else if (strncmp (typeval, "COMP", 4) == 0)
	info.img_types.push_back (IMAGIC_COMP);
      else
	; // throw an error
    }
    fclose (fp);
  }
  return (info);
}

/* Instantiation ----------------------------------------------------------- */
template class ImageImagicT<float>;
template class ImageImagicT<double>;
template class ImageImagicT<complex<double> >;
template bool ImagicWriteImagicFile (const FileName &,
				     const vector<ImageT<float> *>&,
				     ImageImagicType);
template bool ImagicWriteImagicFile (const FileName &,
				     const vector<ImageT<double> *>&,
				     ImageImagicType);
