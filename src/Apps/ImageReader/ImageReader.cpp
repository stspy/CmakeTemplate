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

int main()
{
	string strImgFileName = "E:/Features/test_1.bmp";

	Mat iImgData;

	// get image
	iImgData = imread(strImgFileName);

	if (iImgData.empty())
	{
		return 0;
	}

	imshow("Image", iImgData);
//  imwrite("E:/Features/test_2.bmp", iImgData);

	return 0;
}
