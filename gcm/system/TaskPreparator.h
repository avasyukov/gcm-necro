#ifndef _GCM_TASK_PREPARATOR_H
#define _GCM_TASK_PREPARATOR_H  1

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif
#include <tinyxml.h>

#include <string>
#include <vector>

using std::string;
using std::vector;

#include "../datatypes/MeshOutline.h"
#include "../datatypes/ElasticNode.h"
#include "../meshtypes/TetrMesh_1stOrder.h"
#include "../methods/TetrMethod_Plastic_1stOrder.h"
#include "../rheotypes/VoidRheologyCalculator.h"
#include "../system/LoggerUser.h"
#include "VTKSnapshotWriter.h"
#include "VoidCollisionDetector.h"
#include "BruteforceCollisionDetector.h"
#include "CollisionDetectorForLayers.h"
#include "TetrMeshSet.h"
#include "GCMStresser.h"
#include "GCMException.h"
#include "DataBus.h"
#include "BorderCondition.h"
#include "./forms/StepPulseForm.h"
#include "./areas/BoxArea.h"
#include "../methods/border/ExternalForceCalculator.h"

class TaskPreparator: protected LoggerUser
{
public:
	TaskPreparator();
	~TaskPreparator();
	string* get_task_preparator_type();
	void set_fixed_elastic_rheology(vector<ElasticNode>* nodes, float la, float mu, float rho, float yield_limit);
	void set_fixed_elastic_rheology(vector<ElasticNode>* nodes, MeshOutline* box, float la, float mu, float rho, float yield_limit);
	void set_border_condition(vector<ElasticNode>* nodes, BorderCondition* bc);
	void check_rheology_loaded(vector<ElasticNode>* nodes);

	/**
	 * Do not use it anymore
	 */
	int load_task( string task_file, string zones_file, string data_dir, 
				int* snap_num, int* steps_per_snap);
	/**
	 * Only of of output parameter pSnapStepInterval and pSnapTimeInterval
	 * will be set depending on the configuration. Both of them can't be determined in the same time.
	 * @param task_file
	 * @param zones_file
	 * @param data_dir
	 * @param pTotalStepNum pointer to the filed to receive total the number of steps to take
	 * @param pSnapStepInterval pointer to the filed to receive the number of steps to take before taking next snapshot
	 * @param pSnapTimeInterval pointer to the filed to receive snapshot time interval
	 * @return
	 */
	int load_task (
			string task_file, string zones_file, string data_dir,
			int* pTotalStepNum, int* pSnapStepInterval, float* pSnapTimeInterval);
protected:
	int load_zones_info(string zones_file, vector<int>* zones_info);
	string task_preparator_type;
};

#endif
