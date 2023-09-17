/*=========================================================================

  Program: Remove rectangular regions in pixel data using gdcm

  Copyright (c) 2020 Hauke Bartsch

  =========================================================================*/

#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDefs.h"
#include "gdcmDirectory.h"
#include "gdcmGlobal.h"
#include "gdcmIconImageGenerator.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmReader.h"
#include "gdcmStringFilter.h"
#include "gdcmSystem.h"
#include "gdcmWriter.h"
#include "json.hpp"
#include "optionparser.h"
#include <gdcmUIDGenerator.h>

#include <leptonica/allheaders.h>
#include <tesseract/baseapi.h>

#include <dirent.h>
#include <errno.h>
#include <exception>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <chrono>
#include <filesystem>
#include <map>
#include <pthread.h>
#include <stdio.h>
#include <thread>

struct threadparams
{
  const char **filenames;
  size_t nfiles;
  char *scalarpointer;
  std::string outputdir;
  int thread; // number of the thread
  float confidence;
  bool saveMappings;
  // each thread will store here the study instance uid (original and mapped)
  std::map<std::string, std::string> byThreadStudyInstanceUID;
};

void *ReadFilesThread(void *voidparams)
{
  threadparams *params = static_cast<threadparams *>(voidparams);

  // set params for drawing bounding box
  // Set bounding box parameters
  l_uint32 white_color = 0xFFFFFFFF;

  const size_t nfiles = params->nfiles;
  for (unsigned int file = 0; file < nfiles; ++file)
  {
    const char *filename = params->filenames[file];
    // to get the basename of filename only
    // Find the last '/' character in the filename
    const char *last_slash = strrchr(filename, '/');
    std::string name_without_extension = "test";
    if (last_slash != NULL)
    {
      // Extract the part after the last '/'
      const char *dicom_image_name = last_slash + 1;
      std::string dicom_image_name_str = dicom_image_name; // Convert to string
      size_t dot_position = dicom_image_name_str.find_last_of('.'); // Find the last dot (.) in the string
      if (dot_position != std::string::npos)
      {
          name_without_extension = dicom_image_name_str.substr(0, dot_position); // Extract name without extension
          fprintf(stdout, "File name: %s\n", name_without_extension.c_str()); // Print the name without extension
      }
      else
      {
          fprintf(stdout, "File name: %s\n", dicom_image_name); // No extension found, print the original name
      }
    }

    // std::cerr << filename << std::endl;
    fprintf(stdout, "Start with %s\n", filename);

    gdcm::ImageReader reader;
    // gdcm::Reader reader;
    reader.SetFileName(filename);
    try
    {
      if (!reader.Read())
      {
        std::cerr << "Failed to read: \"" << filename << "\" in thread " << params->thread << std::endl;
        continue;
      }
    }
    catch (...)
    {
      std::cerr << "Failed to read: \"" << filename << "\" in thread " << params->thread << std::endl;
      continue;
    }

    // const gdcm::Image &image = reader.GetImage();
    // if we have the image here we can anonymize now and write again
    gdcm::Anonymizer anon;
    gdcm::File &fileToAnon = reader.GetFile();
    anon.SetFile(fileToAnon);

    gdcm::DataSet &ds = fileToAnon.GetDataSet();

    gdcm::StringFilter sf;
    sf.SetFile(fileToAnon);

    // get some leafs!
    int a = strtol("0020", NULL, 16); // Series Instance UID
    int b = strtol("000E", NULL, 16);
    std::string seriesdirname = sf.ToString(gdcm::Tag(a, b)); // always store by series

    a = strtol("0008", NULL, 16); // SOP Instance UID
    b = strtol("0018", NULL, 16);
    std::string filenamestring = sf.ToString(gdcm::Tag(a, b));

    a = strtol("0020", NULL, 16); // Study Instance UID
    b = strtol("000D", NULL, 16);
    std::string studyinstanceuid = sf.ToString(gdcm::Tag(a, b));

    a = strtol("0008", NULL, 16); // Series Description
    b = strtol("103E", NULL, 16);
    std::string seriesdescription = sf.ToString(gdcm::Tag(a, b));

    a = strtol("0008", NULL, 16); // Study Description
    b = strtol("1030", NULL, 16);
    std::string studydescription = sf.ToString(gdcm::Tag(a, b));

    gdcm::Trace::SetDebug(true);
    gdcm::Trace::SetWarning(true);
    gdcm::Trace::SetError(true);

    int HEIGHT = atoi(sf.ToString(gdcm::Tag(0x0028, 0x0010)).c_str()); // acquisition matrix
    int WIDTH = atoi(sf.ToString(gdcm::Tag(0x0028, 0x0011)).c_str());  // acquisition matrix
    // gdcm::Image *im = new gdcm::Image();
    gdcm::Pixmap &im = reader.GetPixmap(); // is this color or grayscale????
    if (im.AreOverlaysInPixelData())
    { // we can also have curves in here ... what about curves?
      for (int i = 0; i < im.GetNumberOfOverlays(); i++)
      {
        fprintf(stdout, "Warning: we have found an overlay, we will remove it here.\n");
        im.RemoveOverlay(i); // TODO: overlays can hide some text - they are not a good way to anonymize a file
      }
    }
    gdcm::Image gimage = reader.GetImage();

    // -------------------- Create a RGBA image from the DICOM image --------------------------------

    PIX *pixs = pixCreate(WIDTH, HEIGHT, 32); // rgba colors
    size_t length = gimage.GetBufferLength();
    fprintf(stdout, "%ld buffer length size of a single image is: %dx%d\n", length, HEIGHT, WIDTH);
    std::vector<char> vbuffer;
    vbuffer.resize(gimage.GetBufferLength());
    char *buffer = &vbuffer[0];
    if (!gimage.GetBuffer(buffer))
    {
      fprintf(stderr, "Could not get buffer for image data\n");
    }
    // we have a problem with one file:
    bool anyNonZero = false;
    for (int i = 0; i < length; i++)
    {
      if (buffer[i] < buffer[0] || buffer[i] > buffer[0])
      {
        anyNonZero = true;
        break;
      }
    }
    if (anyNonZero)
    {
      fprintf(stdout, "FOUND IMAGE INFORMATION (non-zero char)\n");
    }
    else
    {
      fprintf(stdout, "NO IMAGE INFORMATION FOUND!\n");
    }

    if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::RGB)
    {
      if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
      { 
        fprintf(stdout, "We found RGB data with 8bit! HERE\n");
        unsigned char *ubuffer = (unsigned char *)buffer;
        for (int i = 0; i < HEIGHT; i++)
        {
          for (int j = 0; j < WIDTH; j++)
          {
            l_uint32 val;
            l_int32 red = 0;
            l_int32 green = 0;
            l_int32 blue = 0;
            // if im is grayscale or color
            red = ubuffer[(i * WIDTH + j) * 3 + 0];
            green = ubuffer[(i * WIDTH + j) * 3 + 1];
            blue = ubuffer[(i * WIDTH + j) * 3 + 2];
            composeRGBPixel(red, green, blue, &val);
            pixSetPixel(pixs, j, i, val);
          }
        }
      }
      else
      {
        // color data with not 8bit per channel!
        fprintf(stderr, "Error: found color data with non-8bit per channel! Unsupported!\n");
      }
    }
    else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::MONOCHROME2)
    {
      // we can have 8bit or 16bit grayscales here
      if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
      {
        fprintf(stdout, "We found MONOCHROME2 data with 8bit!\n");
        unsigned char *ubuffer = (unsigned char *)buffer;
        for (int i = 0; i < HEIGHT; i++)
        {
          for (int j = 0; j < WIDTH; j++)
          {
            l_uint32 val;
            l_int32 red = 0;
            l_int32 green = 0;
            l_int32 blue = 0;
            // if im is grayscale or color
            red = ubuffer[i * WIDTH + j];
            green = ubuffer[i * WIDTH + j];
            blue = ubuffer[i * WIDTH + j];
            composeRGBPixel(red, green, blue, &val);
            pixSetPixel(pixs, j, i, val);
          }
        }
      }
      else if (gimage.GetPixelFormat() == gdcm::PixelFormat::INT16)
      { // have not seen an example yet
        short *buffer16 = (short *)buffer;
        fprintf(stdout, "We found MONOCHROME2 data with 16bit!\n");
        for (int i = 0; i < HEIGHT; i++)
        {
          for (int j = 0; j < WIDTH; j++)
          {
            l_uint32 val;
            l_int32 red = 0;
            l_int32 green = 0;
            l_int32 blue = 0;
            // if im is grayscale or color, PixelRepresentation is 1, so we have signed values ->
            // 2complement
            red = (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255);
            green = (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255);
            blue = (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255);
            composeRGBPixel(red, green, blue, &val);
            pixSetPixel(pixs, j, i, val);
          }
        }
      }
      else if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT16)
      { 
        // pixel is 0
        // pixel representation is 0 -> unsigned short
        unsigned short *buffer16 = (unsigned short *)buffer;
        fprintf(stdout, "We found MONOCHROME2 data with 16bit (unsigned short %dx%d)!\n", HEIGHT, WIDTH);
        for (int i = 0; i < HEIGHT; i++)
        {
          for (int j = 0; j < WIDTH; j++)
          {
            l_uint32 val;
            l_int32 red = 0;
            l_int32 green = 0;
            l_int32 blue = 0;
            // if im is grayscale or color
            int v = floor((((double)buffer16[i * WIDTH + j] - gimage.GetPixelFormat().GetMin()) / (float)gimage.GetPixelFormat().GetMax()) * 255);
            red = v; // (unsigned char)std::min(255, (buffer16[i * WIDTH + j]) / 255);
            green = v; // (unsigned char)std::min(255, (buffer16[i * WIDTH + j]) / 255);
            blue = v;  // (unsigned char)std::min(255, (buffer16[i * WIDTH + j]) / 255);
            composeRGBPixel(red, green, blue, &val);
            pixSetPixel(pixs, j, i, val);
          }
        }
      }
      else
      {
        fprintf(stderr, "unknown pixel format in input... nothing is done\n");
      }
    }
    else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::YBR_FULL_422)
    {
      if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
      {
        fprintf(stdout, "Found PhotometricInterpretation YBR_FULL_422 (UINT8)\n");
        unsigned char *ubuffer = (unsigned char *)buffer;
        for (int i = 0; i < HEIGHT; i++)
        {
          for (int j = 0; j < WIDTH; j++)
          {
            l_uint32 val;
            l_int32 red = 0;
            l_int32 green = 0;
            l_int32 blue = 0;
            // if im is grayscale or color
            unsigned char a = ubuffer[(i * WIDTH + j) * 3 + 0];
            unsigned char b = ubuffer[(i * WIDTH + j) * 3 + 1];
            unsigned char c = ubuffer[(i * WIDTH + j) * 3 + 2];

            int R = 38142 * (a - 16) + 52298 * (c - 128);
            int G = 38142 * (a - 16) - 26640 * (c - 128) - 12845 * (b - 128);
            int B = 38142 * (a - 16) + 66093 * (b - 128);

            R = (R + 16384) >> 15;
            G = (G + 16384) >> 15;
            B = (B + 16384) >> 15;

            if (R < 0)
              R = 0;
            if (G < 0)
              G = 0;
            if (B < 0)
              B = 0;
            if (R > 255)
              R = 255;
            if (G > 255)
              G = 255;
            if (B > 255)
              B = 255;
            red = R;
            green = G;
            blue = B;

            composeRGBPixel(red, green, blue, &val);
            pixSetPixel(pixs, j, i, val);
          }
        }
      }
    }
    else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::MONOCHROME1)
    {
      // Check if it's an 8-bit or 16-bit grayscale image
      if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
      {
          fprintf(stdout, "We found MONOCHROME1 data with 8bit!\n");
          unsigned char *ubuffer = (unsigned char *)buffer;
          for (int i = 0; i < HEIGHT; i++)
          {
              for (int j = 0; j < WIDTH; j++)
              {
                  l_uint32 val;
                  l_int32 red = 0;
                  l_int32 green = 0;
                  l_int32 blue = 0;
                  // if im is grayscale or color
                  red = 255 - ubuffer[i * WIDTH + j]; // Invert pixel values
                  green = 255 - ubuffer[i * WIDTH + j]; // Invert pixel values
                  blue = 255 - ubuffer[i * WIDTH + j]; // Invert pixel values
                  composeRGBPixel(red, green, blue, &val);
                  pixSetPixel(pixs, j, i, val);
              }
          }
      }
      else if (gimage.GetPixelFormat() == gdcm::PixelFormat::INT16)
      {
        short *buffer16 = (short *)buffer;
        fprintf(stdout, "We found MONOCHROME1 data with 16bit!\n");
        for (int i = 0; i < HEIGHT; i++)
        {
            for (int j = 0; j < WIDTH; j++)
            {
                l_uint32 val;
                l_int32 red = 0;
                l_int32 green = 0;
                l_int32 blue = 0;
                // if im is grayscale or color, PixelRepresentation is 1, so we have signed values -> 2complement
                red = 255 - (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255); // Invert and scale pixel values
                green = 255 - (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255); // Invert and scale pixel values
                blue = 255 - (unsigned char)std::min(255, (32768 + buffer16[i * WIDTH + j]) / 255); // Invert and scale pixel values
                composeRGBPixel(red, green, blue, &val);
                pixSetPixel(pixs, j, i, val);
            }
        }
      }
      else if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT16)
      {
          unsigned short *buffer16 = (unsigned short *)buffer;
          fprintf(stdout, "We found MONOCHROME1 data with 16bit (unsigned short %dx%d)!\n", HEIGHT, WIDTH);
          for (int i = 0; i < HEIGHT; i++)
          {
              for (int j = 0; j < WIDTH; j++)
              {
                  l_uint32 val;
                  l_int32 red = 0;
                  l_int32 green = 0;
                  l_int32 blue = 0;
                  // Calculate the pixel value within the 0-255 range based on the ratio of the pixel value to the maximum value
                  int v = floor(((double)(buffer16[i * WIDTH + j]) / (float)gimage.GetPixelFormat().GetMax()) * 255);
                  // Invert the pixel value for MONOCHROME1
                  v = 255 - v;
                  red = v;
                  green = v;
                  blue = v;
                  composeRGBPixel(red, green, blue, &val);
                  pixSetPixel(pixs, j, i, val);
              }
          }
      }
    }
    else
    {
      printf("PhotometricInterpretation: %s\n", gimage.GetPhotometricInterpretation().GetString());
      fprintf(stderr, "Error: cannot process this PhotometricInterpretation.\n");
      continue;
    }

    // -------------------- Preprocess the RGBA image before doing OCR thing --------------------------------
    // Invert the colors of the image
    PIX *invertedPix = pixInvert(NULL, pixs);

    // Specify the padding width in pixels
    int paddingWidth = 10; // Adjust this value as needed

    // Load the image with padding
    PIX *pixWithPadding = pixAddBorder(invertedPix, paddingWidth, 0xFFFFFFFF); // White padding

    // -------------------- OCR the RGBA image using tesseract --------------------------------
    std::vector<std::string> safeList = {"Patient", "Name", "Study", "Protocol", "Date", "A", "P", "I", "L", "R", "H"};

    tesseract::TessBaseAPI *api = new tesseract::TessBaseAPI();
    api->Init(NULL, "eng"); // set detect only English
    api->SetPageSegMode(static_cast<tesseract::PageSegMode>(3));

    api->SetImage(pixWithPadding);
    api->SetSourceResolution(70); // set the resolution to 70 dpi

    api->Recognize(0);

    // for debugging write out the pix
    std::string origin_dicom_png = "./data/tmp/" + name_without_extension + "_origin.png";
    std::string thresholded_dicom_png = "./data/tmp/" + name_without_extension + "_thresholded.png";

    if (1)
    {
      pixWrite(origin_dicom_png.c_str(), invertedPix, IFF_PNG);
      Pix *page_pix = api->GetThresholdedImage();
      pixWrite(thresholded_dicom_png.c_str(), page_pix, IFF_PNG);
    }

    tesseract::ResultIterator *ri = api->GetIterator();
    tesseract::PageIteratorLevel level = tesseract::RIL_WORD;
    int counter = 0;
    if (ri != 0)
    {
      do
      {
        if (ri->Empty(level))
        {
          fprintf(stdout, "ignore this level, its empty\n");
          continue;
        }
        const char *word = ri->GetUTF8Text(level);
        const char *word_recognition_language = ri->WordRecognitionLanguage();
        float conf = ri->Confidence(level); // we don't care

        // -------------------------- Get the bounding box and scale that to original image ------------------
        int x1, y1, x2, y2;
        ri->BoundingBox(level, &x1, &y1, &x2, &y2);
        // Scale the bounding box coordinates to the non-padded image
        x1 -= paddingWidth;
        x2 -= paddingWidth;
        y1 -= paddingWidth;
        y2 -= paddingWidth;

        // add some margin and make the selection bigger
        int margin = 0;
        x1 -= margin;
        x2 += margin;
        y1 -= margin;
        y2 += margin;
        if (x1 < 0)
          x1 = 0;
        if (x2 > WIDTH)
          x2 = WIDTH;
        if (y1 < 0)
          y1 = 0;
        if (y2 > HEIGHT)
          y2 = HEIGHT;

        
        // -------------------------- Save the result to a json file ------------------
        if (params->saveMappings)
        {
          // if we store the results we can write them into the thread storage
          char numObjects[11];
          snprintf(numObjects, 11, "%04d", counter++);
          std::string key = filenamestring + "_" + numObjects;
          nlohmann::json info = nlohmann::json::object();
          info["word"] = word ? std::string(word) : std::string("");
          info["confidence"] = conf;
          info["word_recognition_language"] = word_recognition_language ? std::string(word_recognition_language) : std::string("");
          info["word_is_from_dictionary"] = ri->WordIsFromDictionary();
          info["word_is_number"] = ri->WordIsNumeric();
          info["bounding_box"] = nlohmann::json::object({{"x1", x1}, {"y1", y1}, {"x2", x2}, {"y2", y2}});
          info["SOPInstanceUID"] = filenamestring;
          info["SeriesInstanceUID"] = seriesdirname;
          info["StudyInstanceUID"] = studyinstanceuid;
          info["SeriesDescription"] = seriesdescription;
          info["StudyDescription"] = studydescription;
          info["filename"] = filename;
          std::string value = info.dump();

          params->byThreadStudyInstanceUID.insert(std::pair<std::string, std::string>(key, value)); // should only add this pair once
        }

        // -------------------------- Filter the detected word based on detect confidence and some other rules ------------------
        // we can check against a safe list here
        if (std::find(safeList.begin(), safeList.end(), word) != safeList.end())
        {
          printf("skip-word: '%s'; \tconf: %.2f; BoundingBox: %d,%d,%d,%d;\n", word, conf, x1, y1, x2, y2);
          continue; // found a safeList entry, don't do anything
        }

        // check for confidence
        if (conf < params->confidence)
        {
          printf("skip-word - low confidence: '%s'; \tconf: %.2f; BoundingBox: %d,%d,%d,%d;\n", word, conf, x1, y1, x2, y2);
          continue;
        }

        if (strlen(word) == 1)
        {
          printf("skip-word - single character: '%s'; \tconf: %.2f; BoundingBox: %d,%d,%d,%d;\n", word, conf, x1, y1, x2, y2);
          continue;
        }

        // get all valid bounding box to draw on Pix image for demo
        // Iterate over the bounding box and set each pixel to white
        for (int y = y1; y < y2; y++) {
            for (int x = x1; x < x2; x++) {
                pixSetPixel(pixs, x, y, white_color);
            }
        }

        // -------------------- mask out the bounding box (Doing redactor thing) --------------------------------
        printf("word: '%s';  \tconf: %.2f; BoundingBox: %d,%d,%d,%d;\n", word, conf, x1, y1, x2, y2);
        // mask out the bounding box with black
        if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::RGB)
        {
          if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
          { // hopefully always true
            // change values in the buffer
            unsigned char *ubuffer = (unsigned char *)buffer;
            for (int i = y1; i < y2; i++)
            {
              for (int j = x1; j < x2; j++)
              {
                ubuffer[(i * WIDTH + j) * 3 + 0] = 0;
                ubuffer[(i * WIDTH + j) * 3 + 1] = 0;
                ubuffer[(i * WIDTH + j) * 3 + 2] = 0;
              }
            }
          }
        }
        else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::MONOCHROME2)
        {
          // we can have 8bit or 16bit grayscales here
          if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
          {
            unsigned char *ubuffer = (unsigned char *)buffer;
            // fprintf(stdout, "We found MONOCHROME2 data with uint8 8bit!\n");
            for (int i = y1; i < y2; i++)
            {
              for (int j = x1; j < x2; j++)
              {
                ubuffer[i * WIDTH + j] = 0;
              }
            }
          }
          else if (gimage.GetPixelFormat() == gdcm::PixelFormat::INT16)
          { // have not seen an example yet
            short *buffer16 = (short *)buffer;
            // fprintf(stdout, "We found MONOCHROME2 data with 16bit int16!\n");
            for (int i = y1; i < y2; i++)
            {
              for (int j = x1; j < x2; j++)
              {
                buffer16[i * WIDTH + j] = 0;
              }
            }
          }
          else if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT16)
          { 
            unsigned short *buffer16 = (unsigned short *)buffer;
            // fprintf(stdout, "We found MONOCHROME2 data with 16bit (unsigned short)!\n");
            for (int i = y1; i < y2; i++)
            {
              for (int j = x1; j < x2; j++)
              {
                buffer16[i * WIDTH + j] = 0;
              }
            }
          }
        }
        // Son update here for support MONOCHROME1 format
        else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::MONOCHROME1)
        {
            // we can have 8bit or 16bit grayscales here
            if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
            {
                unsigned char *ubuffer = (unsigned char *)buffer;
                fprintf(stdout, "We found MONOCHROME1 data with uint8 8bit!\n");
                for (int i = y1; i < y2; i++)
                {
                    for (int j = x1; j < x2; j++)
                    {
                        ubuffer[i * WIDTH + j] = 255;
                    }
                }
            }
            else if (gimage.GetPixelFormat() == gdcm::PixelFormat::INT16)
            { 
                short *buffer16 = (short *)buffer;
                fprintf(stdout, "We found MONOCHROME1 data with 16bit int16!\n");
                for (int i = y1; i < y2; i++)
                {
                    for (int j = x1; j < x2; j++)
                    {
                        buffer16[i * WIDTH + j] = 32767;
                    }
                }
            }
            else if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT16)
            { 
                unsigned short *buffer16 = (unsigned short *)buffer;
                fprintf(stdout, "We found MONOCHROME1 data with 16bit (unsigned short)!\n");
                for (int i = y1; i < y2; i++)
                {
                    for (int j = x1; j < x2; j++)
                    {
                        buffer16[i * WIDTH + j] = 65535;
                    }
                }
            }
        }
        else if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::YBR_FULL_422)
        {
          if (gimage.GetPixelFormat() == gdcm::PixelFormat::UINT8)
          {                                                   // hopefully always true
                                                              // change values in the buffer
            unsigned char *ubuffer = (unsigned char *)buffer; // but buffer has the wrong encoding now
            for (int i = y1; i < y2; i++)
            {
              for (int j = x1; j < x2; j++)
              { // hope this is ok for YBR data as well
                ubuffer[(i * WIDTH + j) * 3 + 0] = 0;
                ubuffer[(i * WIDTH + j) * 3 + 1] = 0;
                ubuffer[(i * WIDTH + j) * 3 + 2] = 0;
              }
            }
          }
          else
          {
            fprintf(stderr, "ERROR: could not interpret non-UINT8 YBR_FULL_422 data\n");
          }
        }
        else
        {
          fprintf(stdout, "Error: unknown data\n");
        }
        delete[] word;


        // now mask the pixel values
      } while (ri->Next(level));

      // write the redacted image to a file
      std::string redacted_dicom_png = params->outputdir + "/" + name_without_extension + "_redacted.png";
      pixWrite(redacted_dicom_png.c_str(), pixs, IFF_PNG);
    }
    api->End();

    // we need to set the pixel data again that we write, in fileToAnon  (good example
    // https://github.com/malaterre/GDCM/blob/master/Applications/Cxx/gdcmimg.cxx)
    gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
    gdcm::ByteValue *bv = new gdcm::ByteValue();
    bv->SetLength((uint32_t)length); // length here could be strange, something too big for example
    memcpy((void *)bv->GetPointer(), (void *)buffer, length);
    pixeldata.SetValue(*bv);
    if (gimage.GetPhotometricInterpretation() == gdcm::PhotometricInterpretation::YBR_FULL_422)
    { // for YBR_FULL_422
      // we get an error if the transfer syntax is JPEG baseline 1
      gdcm::TransferSyntax ts = gdcm::TransferSyntax::ExplicitVRBigEndian;
      // fileToAnon.GetHeader().SetDataSetTransferSyntax(ts);
      im.SetTransferSyntax(ts);
      // pixeldata.SetVLToUndefined();
    }
    fprintf(stdout, "done...\n");

    // fileToAnon
    im.SetDataElement(pixeldata);
    // reader.SetImage(im);
    // now set the icon as well, we want to rewrite it to make sure we don't have any information in there - overcautions yes!
    if (gimage.GetPhotometricInterpretation() != gdcm::PhotometricInterpretation::RGB)
    {
      gdcm::IconImageGenerator iig;
      iig.AutoPixelMinMax(true);
      iig.SetPixmap(im);
      const unsigned int idims[2] = {64, 64};
      iig.SetOutputDimensions(idims);
      bool b;
      try
      { // this does not catch a float point exception, seems to be causing a crash inside
        b = iig.Generate();
      }
      catch (const std::exception &e)
      {
        fprintf(stderr, "Failed to create icon from image %s\n", e.what());
      }
      if (b)
      {
        const gdcm::IconImage &icon = iig.GetIconImage();
        im.SetIconImage(icon);
      }
    }
    else
    {
      fprintf(stdout, "skip creation of thumb nail image for RGB, can create error on icon generation...\n");
    }

    // ok save the file again
    std::string imageInstanceUID = filenamestring;
    if (imageInstanceUID == "")
    {
      fprintf(stderr, "Error: cannot read image instance uid from %s\n", filename);
      gdcm::UIDGenerator gen;
      imageInstanceUID = gen.Generate();
      filenamestring = imageInstanceUID;
      fprintf(stderr, "Created a random image instance uid: %s\n", imageInstanceUID.c_str());
    }
    std::string fn = params->outputdir + "/" + name_without_extension + "_redacted.dcm";
    // if (1)
    // { // always store results by series directory
    //   // use the series instance uid as a directory name
    //   std::string dn = params->outputdir + "/" + seriesdirname;
    //   struct stat buffer;
    //   if (!(stat(dn.c_str(), &buffer) == 0))
    //   {
    //     mkdir(dn.c_str(), 0777);
    //   } 
    //   fn = params->outputdir + "/" + seriesdirname + "/" + name_without_extension + "_redacted.dcm";
    // }

    fprintf(stdout, "[%d] write to file: %s\n", params->thread, fn.c_str());
    std::string outfilename(fn);

    // save the file again to the output
    gdcm::ImageWriter writer;
    writer.SetFile(fileToAnon);
    writer.SetImage(im);
    writer.SetFileName(outfilename.c_str());
    try
    {
      if (!writer.Write())
      {
        fprintf(stderr, "Error [#file: %d, thread: %d] writing file \"%s\" to \"%s\".\n", file, params->thread, filename, outfilename.c_str());
      }
    }
    catch (const std::exception &ex)
    {
      std::cout << "Caught exception \"" << ex.what() << "\"\n";
    }

    // release resources
    // delete[] buffer;
    // pixDestroy(&pixs);
    // pixDestroy(&invertedPix);
    // pixDestroy(&pixWithPadding);

  }
  return voidparams;
}

void ShowFilenames(const threadparams &params)
{
  std::cout << "start" << std::endl;
  for (unsigned int i = 0; i < params.nfiles; ++i)
  {
    const char *filename = params.filenames[i];
    std::cout << filename << std::endl;
  }
  std::cout << "end" << std::endl;
}

void ReadFiles(size_t nfiles, const char *filenames[], const char *outputdir, int numthreads, float confidence, std::string storeMappingAsJSON)
{
  // \precondition: nfiles > 0
  assert(nfiles > 0);

  // lets change the DICOM dictionary and add some private tags - this is still not sufficient to be
  // able to write the private tags
  gdcm::Global gl;
  if (gl.GetDicts().GetPrivateDict().FindDictEntry(gdcm::Tag(0x0013, 0x0010)))
  {
    gl.GetDicts().GetPrivateDict().RemoveDictEntry(gdcm::Tag(0x0013, 0x0010));
  }
  gl.GetDicts().GetPrivateDict().AddDictEntry(gdcm::Tag(0x0013, 0x0010),
                                              gdcm::DictEntry("Private Creator Group CTP-LIKE", "0x0013, 0x0010", gdcm::VR::LO, gdcm::VM::VM1));

  if (gl.GetDicts().GetPrivateDict().FindDictEntry(gdcm::Tag(0x0013, 0x1010)))
  {
    gl.GetDicts().GetPrivateDict().RemoveDictEntry(gdcm::Tag(0x0013, 0x1010));
  }
  gl.GetDicts().GetPrivateDict().AddDictEntry(gdcm::Tag(0x0013, 0x1010), gdcm::DictEntry("ProjectName", "0x0013, 0x1010", gdcm::VR::LO, gdcm::VM::VM1));

  if (gl.GetDicts().GetPrivateDict().FindDictEntry(gdcm::Tag(0x0013, 0x1013)))
  {
    gl.GetDicts().GetPrivateDict().RemoveDictEntry(gdcm::Tag(0x0013, 0x1013));
  }
  gl.GetDicts().GetPrivateDict().AddDictEntry(gdcm::Tag(0x0013, 0x1013), gdcm::DictEntry("SiteID", "0x0013, 0x1013", gdcm::VR::LO, gdcm::VM::VM1));

  if (gl.GetDicts().GetPrivateDict().FindDictEntry(gdcm::Tag(0x0013, 0x1012)))
  {
    gl.GetDicts().GetPrivateDict().RemoveDictEntry(gdcm::Tag(0x0013, 0x1012));
  }
  gl.GetDicts().GetPrivateDict().AddDictEntry(gdcm::Tag(0x0013, 0x1012), gdcm::DictEntry("SiteName", "0x0013, 0x1012", gdcm::VR::LO, gdcm::VM::VM1));

  if (nfiles <= numthreads)
  {
    numthreads = 1; // fallback if we don't have enough files to process
  }

  const unsigned int nthreads = numthreads; // how many do we want to use?
  threadparams params[nthreads];

  pthread_t *pthread = new pthread_t[nthreads];

  // There is nfiles, and nThreads
  assert(nfiles >= nthreads);
  const size_t partition = nfiles / nthreads;
  for (unsigned int thread = 0; thread < nthreads; ++thread)
  {
    params[thread].filenames = filenames + thread * partition;
    params[thread].outputdir = outputdir;
    params[thread].nfiles = partition;
    params[thread].thread = thread;
    params[thread].confidence = confidence;
    params[thread].saveMappings = false;
    if (storeMappingAsJSON.length() > 0)
    {
      params[thread].saveMappings = true; // store the keys in the params section for later export
    }

    if (thread == nthreads - 1)
    {
      // There is slightly more files to process in this thread:
      params[thread].nfiles += nfiles % nthreads;
    }
    assert(thread * partition < nfiles);
    int res = pthread_create(&pthread[thread], NULL, ReadFilesThread, &params[thread]);
    if (res)
    {
      std::cerr << "Unable to start a new thread, pthread returned: " << res << std::endl;
      assert(0);
    }
  }
  // DEBUG
  size_t total = 0;
  for (unsigned int thread = 0; thread < nthreads; ++thread)
  {
    total += params[thread].nfiles;
  }
  assert(total == nfiles);
  // END DEBUG

  for (unsigned int thread = 0; thread < nthreads; thread++)
  {
    pthread_join(pthread[thread], NULL);
  }

  if (storeMappingAsJSON.length() > 0)
  {
    std::map<std::string, std::string> uidmappings; // make the entries unique by storing them in a map
    for (unsigned int thread = 0; thread < nthreads; thread++)
    {
      for (std::map<std::string, std::string>::iterator it = params[thread].byThreadStudyInstanceUID.begin();
           it != params[thread].byThreadStudyInstanceUID.end(); ++it)
      {
        uidmappings.insert(std::pair<std::string, std::string>(it->first, it->second));
      }
    }
    nlohmann::json ar;
    for (std::map<std::string, std::string>::iterator it = uidmappings.begin(); it != uidmappings.end(); ++it)
    {
      // each value is a json object, to be nice save those as plain json again
      ar[it->first] = nlohmann::json::parse(it->second);
    }

    // if the file exists already, don't store again (maybe its from another thread?)
    std::ofstream jsonfile(storeMappingAsJSON);
    jsonfile << ar;
    jsonfile.flush();
    jsonfile.close();
  }

  delete[] pthread;
}

struct Arg : public option::Arg
{
  static option::ArgStatus Required(const option::Option &option, bool) { return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK; }
  static option::ArgStatus Empty(const option::Option &option, bool) { return (option.arg == 0 || option.arg[0] == 0) ? option::ARG_OK : option::ARG_IGNORE; }
};

enum optionIndex
{
  UNKNOWN,
  HELP,
  INPUT,
  OUTPUT,
  NUMTHREADS,
  CONFIDENCE,
  STOREMAPPING
};
const option::Descriptor usage[] = {{UNKNOWN, 0, "", "", option::Arg::None,
                                     "USAGE: rewritepixel [options]\n\n"
                                     "Options:"},
                                    {HELP, 0, "", "help", Arg::None,
                                     "  --help  \tRewrite DICOM images to remove text. Read DICOM image series and write "
                                     "out an anonymized version of the image data."},
                                    {INPUT, 0, "i", "input", Arg::Required, "  --input, -i  \tInput directory."},
                                    {OUTPUT, 0, "o", "output", Arg::Required, "  --output, -o  \tOutput directory."},
                                    {CONFIDENCE, 0, "c", "confidence", Arg::Required, "  --confidence, -c  \tConfidence threshold (0..100)."},
                                    {NUMTHREADS, 0, "t", "numthreads", Arg::Required, "  --numthreads, -t  \tHow many threads should be used (default 4)."},
                                    {STOREMAPPING, 0, "m", "storemapping", Arg::Required, "  --storemapping, -m  \tStore the detected strings as a JSON file."},
                                    {UNKNOWN, 0, "", "", Arg::None,
                                     "\nExamples:\n"
                                     "  rewritepixel --input directory --output directory\n"
                                     "  rewritepixel --help\n"},
                                    {0, 0, 0, 0, 0, 0}};

// get all files in all sub-directories, does some recursive calls so it might take a while
std::vector<std::string> listFiles(const std::string &path, std::vector<std::string> files)
{
  if (auto dir = opendir(path.c_str()))
  {
    while (auto f = readdir(dir))
    {
      // check for '.' and '..', but allow other names that start with a dot
      if (((strlen(f->d_name) == 1) && (f->d_name[0] == '.')) || ((strlen(f->d_name) == 2) && (f->d_name[0] == '.') && (f->d_name[1] == '.')))
        continue;
      if (f->d_type == DT_DIR)
      {
        std::vector<std::string> ff = listFiles(path + "/" + f->d_name + "/", files);
        // append the returned files to files
        for (int i = 0; i < ff.size(); i++)
        {
          // problem is that we add files here several times if we don't check first.. don't understand why
          if (std::find(files.begin(), files.end(), ff[i]) == files.end())
            files.push_back(ff[i]);
        }
      }

      if (f->d_type == DT_REG)
      {
        // cb(path + f->d_name);
        if (std::find(files.begin(), files.end(), path + "/" + f->d_name) == files.end())
          files.push_back(path + "/" + f->d_name);
      }
    }
    closedir(dir);
  }
  return files;
}

int main(int argc, char *argv[])
{

  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present

  option::Stats stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

  if (parse.error())
    return 1;

  if (options[HELP] || argc == 0)
  {
    option::printUsage(std::cout, usage);
    return 0;
  }

  for (option::Option *opt = options[UNKNOWN]; opt; opt = opt->next())
    std::cout << "Unknown option: " << std::string(opt->name, opt->namelen) << "\n";

  std::string input;
  std::string output;
  int numthreads = 4;
  float confidence = 60.0f; // no confidence is ok
  std::string storeMappingAsJSON = "";
  for (int i = 0; i < parse.optionsCount(); ++i)
  {
    option::Option &opt = buffer[i];
    switch (opt.index())
    {
    case HELP:
      // not possible, because handled further above and exits the program
    case INPUT:
      if (opt.arg)
      {
        fprintf(stdout, "--input '%s'\n", opt.arg);
        input = opt.arg;
      }
      else
      {
        fprintf(stdout, "--input needs a directory specified\n");
        exit(-1);
      }
      break;
    case OUTPUT:
      if (opt.arg)
      {
        fprintf(stdout, "--output '%s'\n", opt.arg);
        output = opt.arg;
      }
      else
      {
        fprintf(stdout, "--output needs a directory specified\n");
        exit(-1);
      }
      break;
    case NUMTHREADS:
      if (opt.arg)
      {
        fprintf(stdout, "--numthreads %d\n", atoi(opt.arg));
        numthreads = atoi(opt.arg);
      }
      else
      {
        fprintf(stdout, "--numthreads needs an integer specified\n");
        exit(-1);
      }
      break;
    case CONFIDENCE:
      if (opt.arg)
      {
        fprintf(stdout, "--confidence %f\n", atof(opt.arg));
        confidence = atoi(opt.arg);
      }
      else
      {
        fprintf(stdout, "--confidence needs an integer specified (0..100)\n");
        exit(-1);
      }
      break;
    case STOREMAPPING:
      if (opt.arg)
      {
        fprintf(stdout, "--storemapping %s\n", opt.arg);
        storeMappingAsJSON = opt.arg;
      }
      else
      {
        fprintf(stdout, "--storemapping needs an path name\n");
        exit(-1);
      }
      break;
    case UNKNOWN:
      // not possible because Arg::Unknown returns ARG_ILLEGAL
      // which aborts the parse with an error
      break;
    }
  }

  // Check if user passed in a single directory - parse all files in all sub-directories
  if (gdcm::System::FileIsDirectory(input.c_str()))
  {
    std::vector<std::string> files;
    files = listFiles(input.c_str(), files);

    const size_t nfiles = files.size();
    const char **filenames = new const char *[nfiles];
    for (unsigned int i = 0; i < nfiles; ++i)
    {
      filenames[i] = files[i].c_str();
      if (std::filesystem::is_symlink(filenames[i]))
      { // replace symlinks with real file names
        filenames[i] = std::filesystem::read_symlink(filenames[i]).c_str();
      }
    }
    if (nfiles == 0)
    {
      fprintf(stdout, "Warning: No files found to process. Skip this step.\n");
      return 0;
    }
    else
    {
      fprintf(stdout, "Info: processing %ld files now.\n", nfiles);
    }
    if (numthreads > nfiles)
      numthreads = nfiles;
    ReadFiles(nfiles, filenames, output.c_str(), numthreads, confidence, storeMappingAsJSON);
    delete[] filenames;
  }
  else
  {
    // its a single file, process that
    const char **filenames = new const char *[1];
    filenames[0] = input.c_str();
    ReadFiles(1, filenames, output.c_str(), 1, confidence, storeMappingAsJSON);
  }

  return 0;
}
