

//#include "stdafx.h"
#include "vtkCamera.h"
#include "vtkGenericRenderWindowInteractor.h"
#include "vtkInteractorStyleJoystickCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLODActor.h"
#include "vtkLight.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPropPicker.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSTLReader.h"
#include "vtkShrinkPolyData.h"

using namespace std;

int main()
{
// 	int nPointNum = 5;
// 
// 	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();
// 	iPoints->InsertNextPoint(0, 0, 0);
// 	iPoints->InsertNextPoint(1, 1, 1);
// 	iPoints->InsertNextPoint(0, 1, 0);
// 	iPoints->InsertNextPoint(0, 1, 1);
// 	iPoints->InsertNextPoint(0, 0, 0);
// 	
// 	vtkSmartPointer<vtkPolyLine> iLine = vtkSmartPointer<vtkPolyLine>::New();
// 	iLine->GetPointIds()->SetNumberOfIds(nPointNum);
// 	for (int i = 0; i < nPointNum; i++)
// 	{
// 		iLine->GetPointIds()->SetId(i, i);
// 	}
// 
// 	vtkSmartPointer<vtkUnstructuredGrid> iGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
// 	iGrid->SetPoints(iPoints);
// 	iGrid->InsertNextCell(iLine->GetCellType(),iLine->GetPointIds());
// 
// 	vtkSmartPointer<vtkDataSetMapper> iMap = vtkSmartPointer<vtkDataSetMapper>::New();
// 	iMap->SetInputData(iGrid);
// 	iMap->SetColorModeToDefault();
// 	iMap->SetScalarModeToUsePointFieldData();
// 	iMap->SelectColorArray("colors");
// 	iMap->SetColorModeToMapScalars();
// 
// 	vtkSmartPointer<vtkActor> iActor = vtkSmartPointer<vtkActor>::New();
// 	iActor->SetMapper(iMap);
// 	iActor->GetProperty()->SetColor(0.8, 0.5, 0.75);
// 
// 	vtkSmartPointer<vtkAssembly> iAssembly = vtkSmartPointer<vtkAssembly>::New();
// 	iAssembly->AddPart(iActor);	
// 
// 	vtkSmartPointer<vtkRenderer> iRen = vtkSmartPointer<vtkRenderer>::New();
// 	iRen->AddActor(iAssembly);
// 	iRen->SetBackground(0, 0, 0);
// 	iRen->ResetCamera();
// 	iRen->GetActiveCamera()->Zoom(1.0);
// 
// 	vtkSmartPointer<vtkRenderWindow> iRenWin = vtkSmartPointer<vtkRenderWindow>::New();
// 	iRenWin->AddRenderer(iRen);
// 	iRenWin->SetSize(500, 500);   
// 
// 	vtkSmartPointer<vtkRenderWindowInteractor> iRenWinIntr = vtkSmartPointer<vtkRenderWindowInteractor>::New();
// 	iRenWinIntr->SetRenderWindow(iRenWin);  
// 	iRenWinIntr->Initialize();
// 	iRenWinIntr->Start();
// 
// 	iRenWin->Render();
// 
// 	return 0;
	vtkRenderer *ren1 = vtkRenderer::New();
	cout<<"aaa"<<endl;
	ren1->GetActiveCamera()->SetClippingRange(0.294421 , 29.4421);
	cout<<"aaa"<<endl;
	ren1->GetActiveCamera()->SetDistance(7.94348);
	cout<<"aaa"<<endl;
	ren1->GetActiveCamera()->SetFocalPoint(-66.9367 , -49.4539 , 258.453);
	ren1->GetActiveCamera()->SetPosition(-67.8091 , -57.3489 , 258.377);
	ren1->GetActiveCamera()->SetViewAngle(20);
	ren1->GetActiveCamera()->SetViewUp(-0.82718 , 0.0860684 , 0.555306);
	ren1->GetActiveCamera()->SetParallelProjection(0);
	ren1->GetActiveCamera()->SetUseHorizontalViewAngle(0);
	ren1->SetBackground(0.1 , 0.2 , 0.4);
	ren1->SetLightFollowCamera(1);
	//创建绘制窗口
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren1);
	renWin->SetSize(1134 , 624);
	//创建交互器
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	iren->SetLightFollowCamera(1);
	//读源对象读取stl数据文件
	vtkSTLReader *part = vtkSTLReader::New();
	part->SetOutput(part->GetOutput());
	part->SetFileName("42400-IDGH.stl");
	//创建过滤器对象，该对象将输入数据集的每个单元向单元质心收缩
	//将会导致相邻单元之间出现裂缝
	vtkShrinkPolyData *shrink = vtkShrinkPolyData::New();
	//将源对象和过滤器连接
	//shrink->SetInput((vtkPolyData *) part->GetOutput());
	shrink->SetInputData((vtkPolyData *) part->GetOutput());
	//设置收缩系数，如果为1，不收缩
	shrink->SetShrinkFactor(0.9);
	//创建映射器对象
	vtkPolyDataMapper *partMapper = vtkPolyDataMapper::New();
	partMapper->SetInputData((vtkPolyData *) shrink->GetOutput());
	partMapper->SetNumberOfPieces(1);
	partMapper->SetScalarRange(0 , 1);
	partMapper->SetColorMode(0);
	partMapper->SetResolveCoincidentTopology(0);
	partMapper->SetScalarMode(0);
	partMapper->SetImmediateModeRendering(0);
	partMapper->SetScalarVisibility(1);
	partMapper->SetUseLookupTableScalarRange(0);
	//创建Props对象(Actor)
	vtkLODActor *partActor = vtkLODActor::New();
	partActor->SetMapper(partMapper);
	partActor->GetProperty()->SetAmbientColor(0.8275 , 0.8275 , 0.8275);
	partActor->GetProperty()->SetColor(0.8275 , 0.8275 , 0.8275);
	partActor->GetProperty()->SetDiffuseColor(0.8275 , 0.8275 , 0.8275);
	partActor->GetProperty()->SetOpacity(1);
	partActor->GetProperty()->SetInterpolation(1);
	partActor->GetProperty()->SetRepresentation(2);
	partActor->GetProperty()->SetBackfaceCulling(0);
	partActor->GetProperty()->SetEdgeVisibility(0);
	partActor->GetProperty()->SetFrontfaceCulling(0);
	partActor->SetOrigin(0 , 0 , 0);
	partActor->SetPosition(0 , 0 , 0);
	partActor->SetScale(1 , 1 , 1);
	partActor->SetVisibility(1);
	//将Actor对象添加到绘制器中
	ren1->AddActor( partActor );
	//绘制
	ren1->ResetCamera();
	ren1->ResetCameraClippingRange();
	renWin->Render();
	iren->Initialize();
	iren->Start();
	//删除对象
	iren->Delete();
	part->Delete();
	partActor->Delete();
	partMapper->Delete();
	ren1->Delete();
	renWin->Delete();
	shrink->Delete();
	return 0;
}

