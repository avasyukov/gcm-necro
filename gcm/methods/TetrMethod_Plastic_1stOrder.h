#ifndef _GCM_TETR_PLASTIC_INTERPOLATION_1ST_ORDER_ROTATE_AXIS_H
#define _GCM_TETR_PLASTIC_INTERPOLATION_1ST_ORDER_ROTATE_AXIS_H  1

#include <gsl/gsl_linalg.h>
#include <vector>
#include <set>
using std::vector;
using std::set;

#include "TetrNumericalMethod.h"
#include "../datatypes/Basis.h"
#include "../datatypes/ElasticNode.h"
#include "../datatypes/ElasticMatrix3D.h"
#include "../system/quick_math.h"
#include "../meshtypes/TetrMesh.h"
#include "../system/TetrMeshSet.h"
#include "../system/SPHConnector.h"
#include "./volume/SimpleVolumeCalculator.h"
#include "./border/FreeBorderCalculator.h"
#include "./border/FixedBorderCalculator.h"
#include "./border/ExternalForceCalculator.h"
#include "./border/ExternalVelocityCalculator.h"
#include "./border/ExternalValuesCalculator.h"
#include "./contact/AdhesionContactCalculator.h"

class GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis : public TetrNumericalMethod
	// TODO may be we should inherit methods from GCM_Tetr_Plastic_Interpolation_1stOrder
{
public:
	GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis();
	~GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis();
	void do_next_part_step(ElasticNode* cur_node, ElasticNode* new_node, float time_step, int stage, TetrMesh* mesh);
	int get_number_of_stages();
	float get_max_lambda(ElasticNode* node, TetrMesh* mesh);
protected:
	// Used for real node
	ElasticMatrix3D* elastic_matrix3d[3];
	// Used for interpolated virtual node in case of contact algorithm
	ElasticMatrix3D* virt_elastic_matrix3d[3];

	SimpleVolumeCalculator* volume_calc;
	FreeBorderCalculator* free_border_calc;
	FixedBorderCalculator* fixed_border_calc;
	ExternalForceCalculator* ext_force_calc;
	ExternalVelocityCalculator* ext_v_calc;
	ExternalValuesCalculator* ext_val_calc;
	AdhesionContactCalculator* adhesion_contact_calc;
	SPHConnector* sph_connector;

	int prepare_node(ElasticNode* cur_node, ElasticMatrix3D* matrixes[], float time_step, int stage, TetrMesh* mesh, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num);
	int prepare_node(ElasticNode* cur_node, ElasticMatrix3D* matrixes[], float time_step, int stage, TetrMesh* mesh, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num, bool debug);

	void prepare_part_step(ElasticNode* cur_node, ElasticMatrix3D* matrix, int stage, int basis_num);
	void drop_deviator(ElasticNode* cur_node, ElasticNode* new_node);
	int create_random_axis(ElasticNode* cur_node, TetrMesh* mesh);

	void log_node_diagnostics(ElasticNode* cur_node, int stage, float outer_normal[], TetrMesh* mesh, int basis_num, ElasticMatrix3D* matrixes[], float time_step, ElasticNode previous_nodes[], int ppoint_num[], bool inner[], float dksi[]);

	// FIXME
	void create_E_matrix(int node_num);

	// These functions are used by create_random_axis() internally
	void create_rotation_matrix(int node_num, float phi, float psi, float teta);
	void create_rotation_matrix_around_normal(int node_num, float x, float y, float z, float teta);
	void create_rotation_matrix_with_normal(int node_num, float x, float y, float z, float teta);
	void find_inverse_matrix(int node_num);

	basis* random_axis;	// New random basises for different nodes
	basis* random_axis_inv;
	int basis_quantity;

	int find_nodes_on_previous_time_layer(ElasticNode* cur_node, int stage, TetrMesh* mesh, float alpha, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num);
	int find_nodes_on_previous_time_layer(ElasticNode* cur_node, int stage, TetrMesh* mesh, float alpha, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num, bool debug);

	quick_math qm_engine;
};

#endif
