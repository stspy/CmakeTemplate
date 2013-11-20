#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/nonfree/nonfree.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace cv;
using namespace std;


const char* keys =
{
	 "{1|  | 0 | camera number}"
};

int main(int argc, const char** argv)
{
	VideoCapture iCap;

	CommandLineParser iParser(argc, argv, keys);
	int nCamNum = iParser.get<int>("1");

	iCap.open(nCamNum);

	if( !iCap.isOpened() )
	{
		cout << "***Could not initialize capturing...***\n";
		cout << "Current parameter's value: \n";
		iParser.printParams();
		return -1;
	}

	Mat iImage, iFrame, iHsvImage, iHueImage;
	bool bPaused = false;

	while (1)
	{
		if( !bPaused )
		{
			iCap >> iFrame;
			if( iFrame.empty() )
				break;
		}

		iFrame.copyTo(iImage);

		if( !bPaused )
		{
			cvtColor(iImage, iHsvImage, CV_BGR2HSV);

			iHueImage.create(iHsvImage.size(), CV_8UC1);
			
			for (int ix = 0; ix < iHsvImage.rows; ix++)
			{
				for (int iy = 0; iy < iHsvImage.cols; iy++)
				{
					iHueImage.at<uchar>(ix, iy) = iHsvImage.at<Vec3b>(ix, iy)[2];
				}
			}

			namedWindow( "Camera", 0 );
			imshow( "Camera", iImage );

			namedWindow( "HSV", 0 );
			imshow( "HSV", iHueImage );

			char c = (char)waitKey(10);
			if( c == 27 )
				break;
		}
	}

	return 0;
}
