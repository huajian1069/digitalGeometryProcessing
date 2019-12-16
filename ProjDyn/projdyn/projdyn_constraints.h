// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// This file contains the Constraint class, which implements ShapeUp
// and Projective Dynamics style constraints.
// First, an abstract Constraint class is defined, which acts as a
// blueprint for specific constraints.
// Then, a few specific constraints are implemented.

#pragma once

#include "projdyn_types.h"
#include "projdyn_common.h"

// Vertices cannot be assigned less weight then this:
constexpr double PROJDYN_MIN_WEIGHT = 1e-6;
// If set to true, the TriangleBendingConstraints prevent the mean curvature
// normal from pointing into the opposite direction as in their original state
constexpr bool PROJDYN_PREVENT_BENDING_FLIPS = false;
// The flipping prevention above is only enforced when the total mean curvature is
// above the value below
constexpr ProjDyn::Scalar PROJDYN_BENDING_FLIP_THRESHOLD = 1e-5;
const float  PI_F=3.14159265358979f;
namespace ProjDyn {
    /* Constraints are grouped into ConstraintGroups, which can have a name, and an additional multiplier.
       This allows for easier accessing from outside. */
    struct ConstraintGroup {
        ConstraintGroup(std::string groupName, const std::vector<ConstraintPtr>& groupConstraints, Scalar groupWeight) {
            name = groupName;
            constraints = groupConstraints;
            weight = groupWeight;
        }

        std::string name = "Unknown";
        std::vector<ConstraintPtr> constraints;
        Scalar weight = 1.;
    };





    /* Blueprint of a Projectiv Dynamics constaint (weight, p, A) */
    class Constraint {
    public:
        /** Constraint constructor.
                vertex_indices - a list of indices of vertices that take part in this constraint
        */
        Constraint(const std::vector<Index>& vertex_indices, Scalar weight) {
            m_vertex_indices = vertex_indices;
            m_weight = (std::max)(PROJDYN_MIN_WEIGHT, weight);
        }

        virtual ~Constraint() {}

        /** Performs the constraint projection
                positions - current positions of the mesh vertices
                projections - a matrix that will be filled by the output of the constraint projection
                            - NOTE: this matrix contains all constraints, and this projection should
                              fill the rows starting from m_constraint_id
        */
        virtual void project(const Positions& positions, Positions& projections) = 0;

        /** Add the constraint to the matrix A that maps vertex positions to the linear part of the constraint projection.
            Specifically, adds a row w_i S_i A_i to the matrix (as triplets), where the notation of the paper is used.
                triplets - A list of all triplets (row, col, entry) that is being built to construct the full matrix A
                currentRow - The first row index of the constraint in the sparse matrix, i.e. the first row in which to add the entries
                    On return this will containt the next free row in the matrix
                sqrtWeight - If true, sqrt(w_i) instead of w_i is used in the formula above
                transpose - If true, the triplets for the transpose of this matrix are returned
        */
        void addConstraint(std::vector<Triplet>& triplets, Index& currentRow, bool sqrtWeight = false, bool transpose = false) {
            m_constraint_id = currentRow;
            std::vector<Triplet> raw_triplets = getTriplets(currentRow);
            for (auto& t : raw_triplets) {
                Scalar v;
                if (sqrtWeight) {
                    v = t.value() * std::sqrt(getWeight());
                }
                else {
                    v = t.value() * getWeight();
                }
                Index row, col;
                if (transpose) {
                    col = t.row();
                    row = t.col();
                }
                else {
                    col = t.col();
                    row = t.row();
                }
                triplets.push_back(Triplet(row, col, v));
            }
            currentRow += getNumConstraintRows();
        }

        /** Returns number of indices of vertices involved in the constraint. */
        Index numIndices() const { return m_vertex_indices.size(); }

        const std::vector<Index>& getIndices() const { return m_vertex_indices; }

        /** Returns the weight of the constraint */
        Scalar getWeight() const { return m_weight * m_weight_mult; }

        /** Change the weight of this constraint. 
            NOTE: this change will only affect the simulation once initializeSystem() is called again. */
        void setWeight(Scalar weight) { m_weight = weight; }

        /** Change the weight multiplier of this constraint.
            NOTE: this change will only affect the simulation once initializeSystem() is called again. */
        void setWeightMultiplier(Scalar mult) { m_weight_mult = mult; }

        virtual ConstraintPtr copy() = 0;

        virtual void set_angle_target(int index, ProjDyn::Scalar angle){;}
    protected:
        /** Add the constraint to the matrix A that maps vertex positions to the linear part of the constraint projection.
            currentRow - the first row index of the constraint in the sparse matrix, i.e. the first row in which to add the entries
            return - a vector of triplets each representing an entry in a sparse matrix
        */
        virtual std::vector<Triplet> getTriplets(Index currentRow) = 0;

        /** Returns the number of rows of the constraint projection auxiliary variable
        */
        virtual Index getNumConstraintRows() = 0;

        /** Indices of the vertices involved in this constraint. */
        std::vector<Index> m_vertex_indices;

        /** Weight for the constraint.*/
        Scalar m_weight = 1.;

        /** Weight multiplier for the constraint. */
        Scalar m_weight_mult = 1.;

        /** Location of this constraint in the linear system. */
        Index m_constraint_id = 0;
    };

    

    class BaseAngleConstraint : public Constraint {
    public:
        BaseAngleConstraint(const std::vector<Index>& edge_vertices,std::vector<int> groupNum, Scalar weight,
            const Positions& positions, bool right_side, 
            std::vector<ProjDyn::Scalar> theta_target, Scalar theta_scope)
            :
            Constraint(edge_vertices, weight)
        {
            // Make sure there are at most two vertices in the edge
            assert(m_vertex_indices.size() == 3);
            ProjDyn::Vector3 left, right, y_normal;
            m_theta_scope = theta_scope;
            m_groupNum = groupNum;
            m_avgPosition.reserve(4);
            updateGroupAvgPosition(positions);

            buildGlobalCoordinate(right_side);

            buildLocalCoordinate(left, right, y_normal);
            m_length.push_back(left.norm());

            extract3DoFanlges(m_yaw_target, m_pitch_target, m_roll_target);
            //m_yaw_target -= 0.5;
            //m_pitch_target -= 0.5;
            //m_roll_target -= 0.5;
            std::cout << " yaw_target: " << m_yaw_target << " pitch_target: " << m_pitch_target << " roll_target: " << m_roll_target << "\n";


            for(int i = 0; i < theta_target.size()-3; i++){
 
                m_theta_target.push_back(extract1DoFanlges(i));
                std::cout << i+4 << "th angle target: " << m_theta_target[i] << std::endl;
                // Compute rest edge length
                m_length.push_back((m_avgPosition[i+2] - m_avgPosition[i+1]).norm());
            }
        }

       

        virtual void project(const Positions& positions, Positions& projection) override {
            // Check for correct size of the projection auxiliary variable;
            assert(projection.rows() > m_constraint_id);
            ProjDyn::Vector3 y_after_2rot, z_after_2rot;
            ProjDyn::Vector3 x_inplane, y_inplane;
            ProjDyn::Vector3 right_final, right;
            Scalar yaw, pitch, roll;
            //build normal vector of first joint plane
            updateGroupAvgPosition(positions);
            extract3DoFanlges(yaw, pitch, roll);
            //std::cout << " yaw: " << yaw << " pitch: " << pitch << " roll: " << roll << "\n";
            //std::cout << " yaw_target: " << m_yaw_target << " pitch_target: " 
            //        << m_pitch_target << " roll_target: " << m_roll_target << "\n";
            //restrict the components in x-y plane
            ProjDyn::Vector3 xy_inplane = ensureAroundTarget(yaw, m_yaw_target, m_x_initial, m_y_initial, "yaw");
            // restrict the components in projection-y plane
            ProjDyn::Vector3 left_final = ensureAroundTarget(pitch, m_pitch_target, xy_inplane, m_z_initial, "pitch");
            projection.row(m_constraint_id) = left_final * m_length[0];


            y_after_2rot = m_z_initial.cross(xy_inplane);
            z_after_2rot = left_final.cross(y_after_2rot);
            ProjDyn::Vector3 y_normal_final = ensureAroundTarget(roll, m_roll_target, y_after_2rot, z_after_2rot, "roll");

            Scalar theta_x = extract1DoFanlges(0);
            //std::cout << "1th angle: " << theta_x << "\n";
            //std::cout << "1th angle target: " << m_theta_target[0] << "\n";
            y_inplane = -left_final.cross(y_normal_final);
            right_final = ensureAroundTarget(theta_x, m_theta_target[0], -left_final, y_inplane, "1th");
            projection.row(m_constraint_id + 1) = right_final * m_length[1];

            x_inplane = -right_final;
            y_inplane = x_inplane.cross(y_normal_final);
            for(int i = 1; i < m_theta_target.size(); i++){
                Scalar theta_x = extract1DoFanlges(i);
           //     std::cout << "2th angle: " << theta_x << "\n";
           //     std::cout << "2th angle target: " << m_theta_target[1] << "\n";
                right_final = ensureAroundTarget(theta_x, m_theta_target[i], x_inplane, y_inplane, "2th");
                projection.row(m_constraint_id + 1 + i) = right_final * m_length[i+1];
                x_inplane = - right_final;
                y_inplane = x_inplane.cross(y_normal_final);
            }

        }

        virtual void set_angle_target(int index, Scalar angle) override{
            assert(index >= 0 && index < 5);
            switch (index){
                case 0:
                    m_yaw_target = angle;
                    break;
                case 1:
                    m_pitch_target = angle;
                    break;
                case 2:
                    m_roll_target = angle;
                    break;
                default:
                    m_theta_target[index-3] = angle;
            }
        }

         void buildGlobalCoordinate(bool right_side){
            // build global coordinate for left and right body
            m_x_initial << -1, 0, 0;
            m_y_initial << 0, 0, 1;
            m_z_initial << 0, 1, 0;
            if(!right_side){
                m_x_initial *= -1;
                m_y_initial *= -1; 
            }
        }

        void buildLocalCoordinate(ProjDyn::Vector3 &left, ProjDyn::Vector3 &right, ProjDyn::Vector3 &y_normal){
            left = m_avgPosition[1] - m_avgPosition[0];
            right = m_avgPosition[2] - m_avgPosition[1];
            y_normal = right.cross(-left).normalized();
        }

        Scalar extract1DoFanlges(int i){
            ProjDyn::Vector3 z_normal, x_inplane, y_inplane, y_normal, left, _;
            buildLocalCoordinate(left, _, y_normal);

            x_inplane = (m_avgPosition[i] - m_avgPosition[i+1]).normalized();
            y_inplane = -y_normal.cross(x_inplane);
            // Compute the normal projection and normalize length to 1
            ProjDyn::Vector3 right = (m_avgPosition[i+2] - m_avgPosition[i+1]);
            ProjDyn::Vector3 projected = (right - right.dot(y_normal) * y_normal).normalized();

            Scalar theta_x = std::acos(right.normalized().dot(x_inplane));
            if(right.dot(y_inplane) < 0)
                theta_x *= -1; 
            return theta_x;
        }

        void extract3DoFanlges(Scalar &yaw, Scalar &pitch, Scalar &roll){
            ProjDyn::Vector3 y_after_2rot, z_after_2rot;
            ProjDyn::Vector3 left, right, y_normal;
            buildLocalCoordinate(left, right, y_normal);
            // yaw
            ProjDyn::Vector3 projected = (left - left.dot(m_z_initial) * m_z_initial).normalized();
            yaw = std::acos(projected.dot(m_x_initial));
            if(projected.dot(m_y_initial) < 0)
                yaw *= -1;

            // pitch
            pitch = PI_F / 2 - std::acos(m_z_initial.dot(left.normalized()));;
            
            // roll
            y_after_2rot = m_z_initial.cross(projected);
            z_after_2rot = left.normalized().cross(y_after_2rot);
            roll = std::acos(y_normal.dot(y_after_2rot));
            if(y_normal.dot(z_after_2rot) < 0)
                roll *= -1;
        }

        ProjDyn::Vector3 ensureAroundTarget(Scalar theta, Scalar target, ProjDyn::Vector3 xAxis, ProjDyn::Vector3 yAxis, std::string str){
            ProjDyn::Vector3 adjusted;
            Scalar theta_scope;
            if(str == "pitch" || str == "roll")
                theta_scope = 0.2;
            else theta_scope = m_theta_scope;
            if(theta < (target - m_theta_scope)){
                std::cout << str + " angle too small " << theta - target  << std::endl;
                adjusted = std::cos(target - theta_scope) * xAxis 
                    + std::sin(target - theta_scope) * yAxis;
            }
            else if(theta > (target + theta_scope)){
                std::cout << str + " angle too big " << theta - target << std::endl;
                adjusted = std::cos(target + theta_scope) * xAxis 
                    + std::sin(target + theta_scope) * yAxis;
            }
            else
                adjusted = std::cos(theta) * xAxis  + std::sin(theta) * yAxis;

            return adjusted;
        }

        void updateGroupAvgPosition(const Positions& positions){
            int num, lower;
            ProjDyn::Vector3 sum;
            for(int i=0; i<4; i++){
                sum << 0.0, 0.0, 0.0;
                if(i == 0) {
                    num = m_groupNum[0];
                    lower = 0;
                }
                else {
                    num = (m_groupNum[i] - m_groupNum[i-1]);
                    lower = m_groupNum[i-1];
                }
                for(int j = lower; j < m_groupNum[i]; j++){
                    sum += positions.row(m_vertex_indices[j]);
                }
                m_avgPosition[i] = sum / num;
            }
        }

        virtual Index getNumConstraintRows() override { return 3; };

        virtual ConstraintPtr copy() {
            return std::make_shared<BaseAngleConstraint>(*this);
        }
    protected:
        Scalar m_theta_scope;
        std::vector<ProjDyn::Scalar> m_theta_target;
        Scalar m_yaw_target, m_pitch_target, m_roll_target;
        std::vector<Scalar> m_length;
        ProjDyn::Vector3 m_x_initial, m_y_initial, m_z_initial;
        std::vector<ProjDyn::Vector3> m_avgPosition;
        std::vector<int> m_groupNum;
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;
            // The current edge is computed as v_k - v_j, where j/k are the
            // indices of the first/second vertex of this edge.
            // In terms of a matrix row A_i and vertex x/y/z-positions v_x/y/z this means
            // A_ij = -1, A_ik = 1 and A_il = 0 for l =/= j,k
            // A_i v_x/y/z = e_x/y/z
            // Thus:
            float num = 1 / m_groupNum[0];
            for(int i = 0; i < m_groupNum[0]; i++){
                triplets.push_back(Triplet(currentRow, m_vertex_indices[i], -num));
            }
            num = 1 / (m_groupNum[1] - m_groupNum[0]);
            for(int i = m_groupNum[0]; i < m_groupNum[1]; i++){
                triplets.push_back(Triplet(currentRow, m_vertex_indices[i], num));
                triplets.push_back(Triplet(currentRow + 1, m_vertex_indices[i], -num));
            }
            num = 1 / (m_groupNum[2] - m_groupNum[1]);
            for(int i = m_groupNum[1]; i < m_groupNum[2]; i++){
                triplets.push_back(Triplet(currentRow + 1, m_vertex_indices[i], num));
                triplets.push_back(Triplet(currentRow + 2, m_vertex_indices[i], -num));
            }
            num = 1 / (m_groupNum[3] - m_groupNum[2]);
            for(int i = m_groupNum[2]; i < m_groupNum[3]; i++){
                triplets.push_back(Triplet(currentRow + 2, m_vertex_indices[i], num));
            }

            return triplets;
        }
    };

    /**	
    As an example for a constraint, the following constraint implements a spring force for an edge of the mesh:
    We have to implement the following methods:
    - the constructor, which initializes the rest data required for this constraint
        (e.g., in this case, the rest length of the spring)
    - the constraint projection, which takes the current positions and computes the "optimal" auxiliary
        projection variable (e.g., in this case, the current edge vector, rescaled to the rest-length)
    - the method that adds the entries for the linear part of the constraint projection
        (e.g., in this case, +1/-1 for the coordinates corresponding to the first/second vertex' coordinate)
    - a method that returns the number of rows for the constraint projection auxiliary variables
        (e.g., in this case, 1, since the auxiliary variable is the edge-vector, with 1 entry per coordinate)
    */
    class EdgeSpringConstraint : public Constraint {
    public:
        EdgeSpringConstraint(const std::vector<Index>& edge_vertices, Scalar weight,
            const Positions& positions)
            :
            Constraint(edge_vertices, weight)
        {
            // Make sure there are at most two vertices in the edge
            assert(m_vertex_indices.size() == 2);
            // Compute rest edge length
            m_rest_length = (positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0])).norm();
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Check for correct size of the projection auxiliary variable;
            assert(projection.rows() > m_constraint_id);
            // Compute the current edge
            projection.row(m_constraint_id) = positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0]);
            // Rescale to rest length
            projection.row(m_constraint_id) /= projection.row(m_constraint_id).norm();
            projection.row(m_constraint_id) *= m_rest_length;
        }

        virtual Index getNumConstraintRows() override { return 1; };

        virtual ConstraintPtr copy() {
            return std::make_shared<EdgeSpringConstraint>(*this);
        }
    protected:
        Scalar m_rest_length = 0;

        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;
            // The current edge is computed as v_k - v_j, where j/k are the
            // indices of the first/second vertex of this edge.
            // In terms of a matrix row A_i and vertex x/y/z-positions v_x/y/z this means
            // A_ij = -1, A_ik = 1 and A_il = 0 for l =/= j,k
            // A_i v_x/y/z = e_x/y/z
            // Thus:
            triplets.push_back(Triplet(currentRow, m_vertex_indices[0], -1));
            triplets.push_back(Triplet(currentRow, m_vertex_indices[1], 1));

            return triplets;
        }
    };

    /**
        Constraint that specifies the position of a group of vertices
    */
    class PositionConstraintGroup : public Constraint {
    public:
        PositionConstraintGroup(std::vector<Index> inds, Scalar weight, const Positions& targetPositions)
            :
            Constraint(inds, weight)
        {
            assert(inds.size() == targetPositions.rows());
            m_target_positions = targetPositions;
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Check for correct size of the projection auxiliary variable;
            assert(projection.rows() >= m_constraint_id + getNumConstraintRows());
            // Set given positions for selected vertices
            for (Index loc_v_ind = 0; loc_v_ind < m_vertex_indices.size(); loc_v_ind++) {
                projection.row(m_constraint_id + loc_v_ind) = m_target_positions.row(loc_v_ind);
            }
        }

        virtual Index getNumConstraintRows() override { return m_vertex_indices.size(); }

        const Positions& getTargetPositions() {
            return m_target_positions;
        }

        void setTargetPositions(const Positions& targetPos) {
            assert(targetPos.rows() == m_target_positions.rows());
            m_target_positions = targetPos;
        }

        virtual ConstraintPtr copy() {
            return std::make_shared<PositionConstraintGroup>(*this);
        }

    protected:
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;
            Index row = 0;
            for (Index v_ind : m_vertex_indices) {
                triplets.push_back(Triplet(currentRow + row, v_ind, 1.));
                row++;
            }
            return triplets;
        }

    private:
        Positions m_target_positions;
    };

    /**
        Constraint that keeps the y position of a vertex above a certain coordinate
    */
    class FloorConstraint : public Constraint {
    public:
        FloorConstraint(Index ind, Scalar weight, Scalar floorHeight, Scalar forceFactor = 1.)
            :
            Constraint({ ind }, weight),
            m_floor_height(floorHeight),
            m_vert_ind(ind),
            m_force_factor(forceFactor)
        {
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Check for correct size of the projection auxiliary variable;
            assert(projection.rows() > m_constraint_id);
            // Set corrected positions for vertices that are below the floor height
            projection.row(m_constraint_id) = positions.row(m_vert_ind);
            if (positions(m_vert_ind, 1) < m_floor_height) {
                projection(m_constraint_id, 1) = (1 + m_force_factor) * m_floor_height - m_force_factor * positions(m_vert_ind, 1);
            }
        }

        virtual Index getNumConstraintRows() override { return 1; }

        virtual ConstraintPtr copy() {
            return std::make_shared<FloorConstraint>(*this);
        }

    protected:
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;
            triplets.push_back(Triplet(currentRow, m_vert_ind, 1.));
            return triplets;
        }

    private:
        Scalar m_floor_height;
        Index m_vert_ind;
        Scalar m_force_factor;
    };

    
    /**
        Tet-strain constraints
    */
    class TetStrainConstraint : public Constraint {
    public:
        TetStrainConstraint(const std::vector<Index>& tetVertices, Scalar weight,
            const Positions& positions, Scalar strain_freedom = 0.)
            :
            Constraint(tetVertices, weight)
        {
            m_strain_freedom = strain_freedom;

            // 3d edges of tet
            Eigen::Matrix<Scalar, 3, 3> edges;
            edges.col(0) = (positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0]));
            edges.col(1) = (positions.row(m_vertex_indices[2]) - positions.row(m_vertex_indices[0]));
            edges.col(2) = (positions.row(m_vertex_indices[3]) - positions.row(m_vertex_indices[0]));

            // Inverse of edges matrix for computation of the deformation gradient
            m_rest_edges_inv = edges.inverse();
            bool didCorrect = false;
            while (!m_rest_edges_inv.allFinite()) {
                std::cout << "Illegal edges in mesh!" << std::endl;
                for (int c = 0; c < 3; c++) {
                    if (edges.col(c).norm() < 1e-12) {
                        for (int r = 0; r < 3; r++) {
                            edges(r, c) = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 1e-11;
                        }
                    }
                }
                for (int r = 0; r < 3; r++) {
                    if (edges.row(r).norm() < 1e-12) {
                        for (int c = 0; c < 3; c++) {
                            edges(r, c) = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 1e-11;
                        }
                    }
                }
                m_rest_edges_inv = edges.inverse();
                didCorrect = true;
            }
            if (didCorrect) std::cout << "Fixed!" << std::endl;
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Compute deformation gradient, clamp its singular values and output
            // corrected deformation gradient as the projection

            // 3d edges of tet
            Eigen::Matrix<Scalar, 3, 3> edges;
            edges.col(0) = (positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0]));
            edges.col(1) = (positions.row(m_vertex_indices[2]) - positions.row(m_vertex_indices[0]));
            edges.col(2) = (positions.row(m_vertex_indices[3]) - positions.row(m_vertex_indices[0]));

            // Compute the deformation gradient (current edges times inverse of original edges)
            Eigen::Matrix<Scalar, 3, 3> F = edges * m_rest_edges_inv;
            // Compute SVD
            Eigen::JacobiSVD<Eigen::Matrix<Scalar, 3, 3>> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Vector3 S = svd.singularValues();
            Scalar detU = svd.matrixU().determinant();
            Scalar detV = svd.matrixV().determinant();

            // Clamp singular values
            S(0) = clamp(S(0), 1. - m_strain_freedom, 1. + m_strain_freedom);
            S(1) = clamp(S(1), 1. - m_strain_freedom, 1. + m_strain_freedom);
            S(2) = clamp(S(2), 1. - m_strain_freedom, 1. + m_strain_freedom);
            // Reverse largest singular value if tet is inverted:
            if (detU * detV < 0.0f) S(2) = -S(2);

            // Compute optimal deformation gradient
            F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

            if (F.hasNaN()) {
                F.setIdentity();
            }

            // Output corrected deformation gradient
            projection.block<3, 3>(m_constraint_id, 0) = F.transpose();
        }

        virtual Index getNumConstraintRows() override { return 3; }

        virtual ConstraintPtr copy() {
            return std::make_shared<TetStrainConstraint>(*this);
        }

    protected:
        Eigen::Matrix<Scalar, 3, 3> m_rest_edges_inv;
        Scalar m_strain_freedom;

        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            // The selection matrix computes the current deformation gradient w.r.t the current positions
            // (i.e. multiplication of the current edges with the inverse of the original edge matrix)
            std::vector<Triplet> triplets;

            for (int coord3d = 0; coord3d < 3; coord3d++) {
                triplets.push_back(Triplet(currentRow + coord3d, m_vertex_indices[0], -(m_rest_edges_inv(0, coord3d) + m_rest_edges_inv(1, coord3d) + m_rest_edges_inv(2, coord3d))));
                triplets.push_back(Triplet(currentRow + coord3d, m_vertex_indices[1], m_rest_edges_inv(0, coord3d)));
                triplets.push_back(Triplet(currentRow + coord3d, m_vertex_indices[2], m_rest_edges_inv(1, coord3d)));
                triplets.push_back(Triplet(currentRow + coord3d, m_vertex_indices[3], m_rest_edges_inv(2, coord3d)));
            }

            return triplets;
        }
    };

    /**
        Triangle-strain constraints
    */
    class TriangleStrainConstraint : public Constraint {
    public:
        TriangleStrainConstraint(const std::vector<Index>& triVertices, Scalar weight,
            const Positions& positions, Scalar strain_freedom = 0.)
            :
            Constraint(triVertices, weight)
        {
            m_strain_freedom = 0.;

            Eigen::Matrix<Scalar, 3, 2> edges, P;
            // 3d edges of triangle
            edges.col(0) = (positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0]));
            edges.col(1) = (positions.row(m_vertex_indices[2]) - positions.row(m_vertex_indices[0]));

            // Projection that embeds these edges isometrically in 2d, in a way that the first edge is aligned to the x-axis
            P.col(0) = edges.col(0).normalized();
            P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();

            // Compute the 2d rest edges
            Eigen::Matrix<Scalar, 2, 2> restEdges = P.transpose() * edges;

            // ... and their inverse
            m_rest_edges_inv = restEdges.inverse();
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Project the current edges isometrically into 2d, compute the deformation gradient there
            // then perform the SVD on the def.grad., clamp the singular values, and project the deformation
            // gradient back.
            Eigen::Matrix<Scalar, 3, 2> edges, P;

            // 3d edges of triangle
            edges.col(0) = (positions.row(m_vertex_indices[1]) - positions.row(m_vertex_indices[0]));
            edges.col(1) = (positions.row(m_vertex_indices[2]) - positions.row(m_vertex_indices[0]));

            // Projection that embeds these edges isometrically in 2d
            P.col(0) = edges.col(0).normalized();
            P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
            // Compute the deformation gradient (current 2d edges times inverse of original 2d edges)
            Eigen::Matrix<Scalar, 2, 2> F = P.transpose() * edges * m_rest_edges_inv;
            // Compute SVD
            Eigen::JacobiSVD<Eigen::Matrix<Scalar, 2, 2>> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Vector2 S = svd.singularValues();
            // Clamp singular values
            S(0) = clamp(S(0), 1. - m_strain_freedom, 1. + m_strain_freedom);
            S(1) = clamp(S(1), 1. - m_strain_freedom, 1. + m_strain_freedom);
            // Compute clamped deformation gradient
            F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

            // Project deformation gradient to 3d
            projection.block<2, 3>(m_constraint_id, 0) = (P * F).transpose();
        }

        virtual Index getNumConstraintRows() override { return 2; }

        virtual ConstraintPtr copy() {
            return std::make_shared<TriangleStrainConstraint>(*this);
        }
    protected:
        Eigen::Matrix<Scalar, 2, 2> m_rest_edges_inv;
        Scalar m_strain_freedom;
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;

            // The selection matrix computes, from the current vertex positions,
            // the deformation gradient (projected to 3d).
            for (int coord2d = 0; coord2d < 2; coord2d++) {
                triplets.push_back(Triplet(currentRow + coord2d, m_vertex_indices[0], -(m_rest_edges_inv(0, coord2d) + m_rest_edges_inv(1, coord2d))));
                triplets.push_back(Triplet(currentRow + coord2d, m_vertex_indices[1], m_rest_edges_inv(0, coord2d)));
                triplets.push_back(Triplet(currentRow + coord2d, m_vertex_indices[2], m_rest_edges_inv(1, coord2d)));
            }

            return triplets;
        }
    };





    /**
        Bending constraints, to be implemented in exercises.cpp
    */
    class BendingConstraint : public Constraint {
    public:
        BendingConstraint(VertexStar vertexStar, Scalar weight, Scalar voronoiArea,
            const Positions& positions, const Triangles& triangles) 
            :
            Constraint(getStarIndices(vertexStar), weight)
        {
            m_vertex_star = vertexStar;

            // Compute cotan weights and collect list of triangles in 1-ring
            Vector cotanWeights(m_vertex_star.size());
            int nb = 0, nb2 = 0;
            std::vector<Index> trisSeen;
            trisSeen.reserve(m_vertex_star.size());
            m_triangles.resize(m_vertex_star.size(), 3);
            for (Edge& e : m_vertex_star) {
                Vector3 edgeT11 = (positions.row(m_vertex_indices[0]) - positions.row(e.vOtherT1));
                Vector3 edgeT12 = (positions.row(e.v2) - positions.row(e.vOtherT1));
                Scalar angle1 = std::acos(edgeT11.dot(edgeT12) / (edgeT11.norm() * edgeT12.norm()));
                Scalar cotWeight = (.5 / std::tan(angle1));
                if (e.t2 > 0) {
                    Vector3 edgeT21 = (positions.row(m_vertex_indices[0]) - positions.row(e.vOtherT2));
                    Vector3 edgeT22 = (positions.row(e.v2) - positions.row(e.vOtherT2));
                    Scalar angle2 = std::acos(edgeT21.dot(edgeT22) / (edgeT21.norm() * edgeT22.norm()));
                    cotWeight += (.5 / std::tan(angle2));
                }

                cotanWeights(nb) = cotWeight / voronoiArea;

                if (std::find(trisSeen.begin(), trisSeen.end(), e.t1) == trisSeen.end()) {
                    trisSeen.push_back(e.t1);
                    m_triangles.row(nb2) = triangles.row(e.t1);
                    nb2++;
                }
                if (e.t2 >= 0 && std::find(trisSeen.begin(), trisSeen.end(), e.t2) == trisSeen.end()) {
                    trisSeen.push_back(e.t2);
                    m_triangles.row(nb2) = triangles.row(e.t2);
                    nb2++;
                }
                nb++;
            }
            if (nb2 < nb) {
                Triangles triTemp = m_triangles.block(0, 0, nb2, 3);
                m_triangles = triTemp;
            }

            // Store cotan weights
            m_cotan_weights = cotanWeights;

            // Compute and store mean-curvature (vector)
            Eigen::Matrix<Scalar, 1, 3> meanCurvatureVector;
            meanCurvatureVector.setZero();
            nb = 0;
            for (Edge e : vertexStar) {
                meanCurvatureVector += (positions.row(m_vertex_indices[0]) - positions.row(e.v2)) * cotanWeights(nb);
                nb++;
            }
            m_rest_mean_curv_vec = meanCurvatureVector;
            m_rest_mean_curv = meanCurvatureVector.norm();

            // Compute and store the dot product of the mean curvature vector with the
            // normal.
            Vector3 triNormal = getTriangleNormal(m_triangles, positions);
            m_dot_with_normal = triNormal.dot(meanCurvatureVector);

            if (m_cotan_weights.hasNaN()) {
                std::cout << "NaN value in cotan weights!" << std::endl;
            }
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            if (m_rest_mean_curv < 1e-12) {
                projection.row(m_constraint_id).setZero();
            }

            // Compute mean curvature vector
            Eigen::Matrix<Scalar, 1, 3> meanCurvatureVector;
            meanCurvatureVector.setZero();
            int nb = 0;
            for (Edge e : m_vertex_star) {
                meanCurvatureVector += (positions.row(m_vertex_indices[0]) - positions.row(e.v2)) * m_cotan_weights(nb);
                nb++;
            }

            // Rescale mean curvature vector to rest-length
            Scalar norm = meanCurvatureVector.norm();
            if (norm < 1e-10) {
                meanCurvatureVector = getTriangleNormal(m_triangles, positions) * m_rest_mean_curv;
            }
            else {
                meanCurvatureVector *= m_rest_mean_curv / norm;
            }

            // If bending flips should be prevented, check if the sign of the dot product
            // with the normal is the same as in the rest configuration, otherwise flip
            // the mean curvature vector.
            if (PROJDYN_PREVENT_BENDING_FLIPS) {
                Vector3 triNormal = getTriangleNormal(m_triangles, positions);
                Scalar dot = triNormal.dot(meanCurvatureVector);
                if (norm > PROJDYN_BENDING_FLIP_THRESHOLD&& dot* m_dot_with_normal < 0) meanCurvatureVector *= -1;
            }

            projection.row(m_constraint_id) = meanCurvatureVector;
        }

        virtual Index getNumConstraintRows() override { return 1; }

        virtual ConstraintPtr copy() {
            return std::make_shared<BendingConstraint>(*this);
        }
    protected:
        VertexStar m_vertex_star;
        Vector m_cotan_weights;
        Scalar m_dot_with_normal;
        Scalar m_rest_mean_curv;
        Eigen::Matrix<Scalar, 1, 3> m_rest_mean_curv_vec;
        Triangles m_triangles;
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            // For bending constraints, the selection matrix computes the current mean curvature
            // vector of each vertex.
            std::vector<Triplet> triplets;
            triplets.push_back(Triplet(currentRow, m_vertex_indices[0], m_cotan_weights.sum()));
            Index nb = 0;
            for (Edge e : m_vertex_star) {
                triplets.push_back(Triplet(currentRow, e.v2, -m_cotan_weights(nb)));
                nb++;
            }
            return triplets;
        }
    };





    class FlatteningConstraint : public BendingConstraint {
    public:
        FlatteningConstraint(VertexStar vertexStar, Scalar weight, Scalar voronoiArea,
            const Positions& positions, const Triangles& triangles)
            : BendingConstraint(vertexStar, weight, voronoiArea, positions, triangles) {};

        virtual void project(const Positions& positions, Positions& projection) override {
            projection.row(m_constraint_id).setZero();
        }

        virtual ConstraintPtr copy() {
            return std::make_shared<FlatteningConstraint>(*this);
        }
    };





    class CenterConstraint : public Constraint {
    public:
        CenterConstraint(const std::vector<Index>& vertex_indices, Scalar weight, Vector3 center)
            :
            Constraint(vertex_indices, weight)
        {
            m_center = center;
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            projection.row(m_constraint_id) = m_center.transpose();
        }

        virtual Index getNumConstraintRows() override { return 1; }

        void setCenter(const Vector3& center) {
            m_center = center;
        }

        const Vector3& getCenter() {
            return m_center;
        }

        virtual ConstraintPtr copy() {
            return std::make_shared<CenterConstraint>(*this);
        }
    protected:
        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            // Compute the average of all positions listed in m_vertex_indices
            std::vector<Triplet> triplets;
            Scalar fac = m_vertex_indices.size();
            for (Index ind : m_vertex_indices) {
                triplets.push_back(Triplet(currentRow, ind, 1/fac));
            }
            return triplets;
        }

        Vector3 m_center;
    };




    class IsometricOneRing : public Constraint {
    public:
        IsometricOneRing(const std::vector<Index>& ring_vertices, Scalar weight,
            const Positions& positions)
            :
            Constraint(ring_vertices, weight)
        {
            // Initialize 1-Ring constraint
            // [Add code here!]
        }

        virtual void project(const Positions& positions, Positions& projection) override {
            // Compute best fit of original one ring edges to current one ring edges
            // Fill these edges into the projection matrix, starting from the row
            // given by m_constraint_id
            // [Add code here!]
        }

        virtual Index getNumConstraintRows() override { 
            return 0; // [Return correct value here!]; 
        };

        virtual ConstraintPtr copy() {
            return std::make_shared<IsometricOneRing>(*this);
        }
    protected:
        // [Add member variables that you need for this constraint here]

        virtual std::vector<Triplet> getTriplets(Index currentRow) override {
            std::vector<Triplet> triplets;
            // Compute the k edges of the current 1-Ring
            // [Add code here!]
            return triplets;
        }
    };
    
}