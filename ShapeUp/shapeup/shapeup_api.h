// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// This file contains an API that bridges the gap between simulator
// and the Viewer to enable ShapeUp-style shape deformation.
// It provides GUI elements to change and add constraints, reset the
// mesh or export it, and more.
// Futhermore it provides call-backs for mouse-events from
// the viewer that handle the click + drag behavior on the mesh.

#pragma once

#include "projdyn.h"
#include "viewer.h"
#include "projdyn_api.h"
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/button.h>


/* begin of modification */
#include <fstream>
#define THETA_SCOPE 0.1         // tolerance scope: (target angle +/- THETA_SCOPE)
#define THETA_TARGET 1.5        // never used
/* end of modification */

using namespace ProjDyn;

// The default weight given to "movable" vertex groups
constexpr ProjDyn::Scalar SHAPEUP_MOVABLE_WEIGHT = 1.;
// The default number of local-global iterations
constexpr ProjDyn::Index SHAPEUP_DEFAULT_NUM_ITS = 10;

class ShapeUpAPI {
public:
    ShapeUpAPI(Viewer* viewer)
        :
        m_pdAPI(viewer)
    {
        // The ProjDynAPI is another layer of API that handles more
        // general communication with the simulator.
        m_pdAPI.setDynamic(false);

        // The viewer displays the mesh and handles mouse events.
        // It also provides a few GUI elements, such as loading meshes
        // and display options.
        m_viewer = viewer;
    }
    
    // Called every time the mesh in the viewer changed.
    // Pass the viewer mesh to the simulator and reset constraints etc.
    bool setMesh() {
        // Reset position constraint groups
        m_movableGroups.clear();

        // Pass the mesh from the viewer to the projective dynamics framework
        if (!m_pdAPI.setMesh()) {
            return false;
        }

        // Store initial vertex positions
        m_initialPositions = m_pdAPI.getPositions();

        // Add default positions constraints
        addCenterConstraint();

        // Reset status of the shape constraints radio button
        if (m_buttonEdgeSprings) {
            if (m_buttonEdgeSprings->pushed()) {
                // Add edge spring constraints
                m_pdAPI.addEdgeSpringConstraints(10.);
            }
            else {
                // Add curvature and area constraints
                m_pdAPI.addBendingConstraints(10.);
                m_pdAPI.addTriangleStrainConstraints(10.);
            }
        }
        else {
            // Add edge spring constraints
            m_pdAPI.addEdgeSpringConstraints(10.);
        }

        // Perform a first local-global solve (and initialize system)
        m_pdAPI.update(true);

        m_isInitialized = true;
        return true;
    }

    // Here, the ShapeUp GUI is created, i.e. buttons to add constraints,
 // reset the shape, export it, etc...
    void initShapeUpGUI() {
        Window* shapeup_win = new Window(m_viewer, "ShapeUp");
        shapeup_win->setPosition(Vector2i(15, 230));
        shapeup_win->setLayout(new GroupLayout());

        Button* b = new Button(shapeup_win, "Reset (restore orig)");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            setMesh();
            updateVertexStatus();
        });

        b = new Button(shapeup_win, "Reset (current -> rest)");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            makeCurrentShapeToRestShape();
            updateVertexStatus();
        });

        b = new Button(shapeup_win, "Export to OBJ");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            // Copy mesh and update positions
            Surface_mesh copy_m(*(m_viewer->getMesh()));
            const Positions& cur_pos = m_pdAPI.getPositions();
            for (const auto& v : copy_m.vertices()) {
                Index ind = v.idx();
                copy_m.position(v).x = cur_pos(ind, 0);
                copy_m.position(v).y = cur_pos(ind, 1);
                copy_m.position(v).z = cur_pos(ind, 2);
            }
            copy_m.write("../export.obj");
        });

        new Label(shapeup_win, "Shape Constraint Type", "sans-bold");
        /*
        b = new Button(shapeup_win, "Edge Springs");
        b->setFlags(Button::RadioButton);
        b->setPushed(false);
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            // Remove all constraints
            m_pdAPI.clearConstraints();
            // Add edge spring constraints
            m_pdAPI.addEdgeSpringConstraints(10.);
            if (m_movableGroups.empty()) {
                // Add center constraint (in absence of position constraints)
                addCenterConstraint();
            }
            else {
                // Add back previous position constraints
                for (auto conGroup : m_movableGroups)
                    m_pdAPI.addConstraints(conGroup);
            }
            updateVertexStatus();
            m_pdAPI.update(true);
            m_usedShapeConstraint = 0;
        });
        m_buttonEdgeSprings = b;
        */
        b = new Button(shapeup_win, "Triangle Area + Curvature");
        b->setFlags(Button::RadioButton);
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            // Remove all constraints
            m_pdAPI.clearConstraints();
            // Add curvature and area constraints
            //m_pdAPI.addBendingConstraints(10.);
            m_pdAPI.addTriangleStrainConstraints(10.);
            if (m_movableGroups.empty()) {
                // Add center constraint (in absence of position constraints)
                addCenterConstraint();
            }
            else {
                // Add back previous position constraints
                for (auto conGroup : m_movableGroups)
                    m_pdAPI.addConstraints(conGroup);
            }
            updateVertexStatus();
            m_pdAPI.update(true);
            m_usedShapeConstraint = 1;
        });

        // Example for an additional constraint that is added to the selected vertices
        b = new Button(shapeup_win, "Iso 1-Rings (EXERCISE)");
        b->setFlags(Button::RadioButton);
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            // Remove all constraints
            m_pdAPI.clearConstraints();
            // Add isometric 1-ring constraints
            addIsometricOneRingConstraints(100.);
            if (m_movableGroups.empty()) {
                // Add center constraint (in absence of position constraints)
                addCenterConstraint();
            }
            else {
                // Add back previous position constraints
                for (auto conGroup : m_movableGroups)
                    m_pdAPI.addConstraints(conGroup);
            }
            updateVertexStatus();
            m_pdAPI.update(true);
            m_usedShapeConstraint = 2;
        });

        new Label(shapeup_win, "Position Constraints", "sans-bold");

        b = new Button(shapeup_win, "Clear position constraints");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            // Remove all constraints
            m_pdAPI.clearConstraints();
            m_movableGroups.clear();
            // Add back shape constraints
            if (m_usedShapeConstraint == 0) {
                m_pdAPI.addEdgeSpringConstraints(10.);
            } else if (m_usedShapeConstraint == 1) {
                m_pdAPI.addBendingConstraints(10.);
                m_pdAPI.addTriangleStrainConstraints(10.);
            } else {
                addIsometricOneRingConstraints(100);
            }            
            // We still need to add a single position constraint
            // otherwise the mesh can fly away
            addCenterConstraint();
            updateVertexStatus();
            m_viewer->clearSelection();
            m_pdAPI.update(true);
        });

        /* begin of modification */
        b = new Button(shapeup_win, "Fix Selected Points");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            const auto& selVerts = m_viewer->getSelectedVertices();
            if (selVerts.size() == 0) return;
            addPositionConstraintGroup(selVerts);
            m_viewer->clearSelection();
            m_pdAPI.update(true);
        });

        b = new Button(shapeup_win, "Cache Selected Group of Vertices");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            const std::vector<Index>& selVerts = m_viewer->getSelectedVertices();
            if(m_avgPosition.size() == 4){
                m_avgPosition.clear();
                m_sortedIndex.clear();
                m_original_edge_length.clear();
                m_groupNum.clear();
            }
            ProjDyn::Vector3 avg;
            avg << 0.0, 0.0, 0.0;
            for(Index vertices : selVerts){
                m_sortedIndex.push_back(vertices);
                avg += m_pdAPI.getPositions().row(vertices);
            }
            avg /= selVerts.size();
            if(m_avgPosition.size() > 0){
                double length = (avg - m_avgPosition.back()).norm();
                m_original_edge_length.push_back(length);
            }
            m_avgPosition.push_back(avg);
            m_groupNum.push_back(selVerts.size());
            if(m_avgPosition.size() == 4){
                for(int i=1; i<m_groupNum.size(); i++)
                    m_groupNum[i] += m_groupNum[i-1]; 
            }
            std::cout << m_groupNum.size() <<"th group of "<< selVerts.size() << " points added" << std::endl;
        });

        b = new Button(shapeup_win, "Display Info of Cached Points");
        b->setCallback([this]() {
            m_avgPosition.reserve(4);
            updateGroupAvgPosition();
            for(int i=0; i<4; i++){
                std::cout << m_avgPosition[i].transpose() << std::endl;
                std::cout << i << "th index " << m_sortedIndex[i] << std::endl;
            }
            std::cout << "theta_scope: " << THETA_SCOPE << " theta_target: " << THETA_TARGET << std::endl;
            ProjDyn::Vector3 m_x_initial, m_y_initial, m_z_initial;
            ProjDyn::Vector3 y_after_2rot, z_after_2rot;
            ProjDyn::Scalar yaw, pitch, roll;
            m_x_initial << -1, 0, 0;
            m_y_initial << 0, 0, 1;
            m_z_initial << 0, 1, 0;

            ProjDyn::Vector3 left = m_avgPosition[1] - m_avgPosition[0];
            ProjDyn::Vector3 right = m_avgPosition[2] - m_avgPosition[1];
            ProjDyn::Vector3 y_normal = right.cross(-left).normalized();
            if(right.dot(-left) < 0)
                y_normal *= -1;

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

            std::cout << " yaw: " << yaw << " pitch: " << pitch << " roll: " << roll << "\n";
            if(m_avgPosition.size() > 1)
                std::cout << "0)length_ratio: current/original = " << (m_avgPosition[1] - 
                    m_avgPosition[0]).norm() / m_original_edge_length[0] << std::endl;      
            for(int i=1; i<m_avgPosition.size()-1; i++){
                double length = (m_avgPosition[i+1] - m_avgPosition[i]).norm();
                double cos = (m_avgPosition[i+1] - m_avgPosition[i]).normalized()
                            .dot((m_avgPosition[i-1] - m_avgPosition[i]).normalized());
                std::cout << i << ")length_ratio: current/original = " <<  length / m_original_edge_length[i]; 
                std::cout << " current angle:  "<< std::acos(cos) << std::endl;
            }
        });
        /* end of modification */

        b = new Button(shapeup_win, "Select Boundary");
        b->setCallback([this]() {
            std::vector<Index> boundInds;
            Surface_mesh* smesh = m_viewer->getMesh();
            for (auto v : smesh->vertices()) {
                if (smesh->is_boundary(v)) boundInds.push_back(v.idx());
            }
            m_viewer->setSelectedVertices(boundInds);
        });

        m_pdAPI.setNumIterations(SHAPEUP_DEFAULT_NUM_ITS);
        Label* iterations_label = new Label(shapeup_win, "Num Loc-Glob Its: ");
        IntBox<int>* iterations_box = new IntBox<int>(shapeup_win, m_pdAPI.getNumIterations());
        iterations_box->setEditable(true);
        iterations_box->setCallback([this](int num_its) {
            if (!m_isInitialized) return;
            m_pdAPI.setNumIterations(num_its);
        });

        new Label(shapeup_win, "Additional constraints on selection:", "sans-bold");

        Widget* panel = new Widget(shapeup_win);
        panel->setLayout(new BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 20));

        // Example for an additional constraint that is added to the selected vertices
        b = new Button(shapeup_win, "Add Flat-Constr. to selection");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            addFlatnessConstraintsToSelection(1.);
            m_pdAPI.update(true);
            m_viewer->clearSelection();
            m_usedShapeConstraint = false;
        });

        /* begin of modification */
        b = new Button(shapeup_win, "Add edge constraints to selected");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            addEdgeSpringConstraintsSelection(1.);
            m_pdAPI.update(true);
            m_viewer->clearSelection();
            m_usedShapeConstraint = false;
        });

        // add constraint to one leg on right side
        b = new Button(shapeup_win, "Add Angle constraints to selected(right) leg");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            addBaseAngleConstraintsSelection(1., true);
            m_pdAPI.update(true);
            m_viewer->clearSelection();
            m_usedShapeConstraint = false;
        });
        // add constraint to one leg on left side
        b = new Button(shapeup_win, "Add Angle constraints to selected(left) leg");
        b->setCallback([this]() {
            if (!m_isInitialized) return;
            addBaseAngleConstraintsSelection(1., false);
            m_pdAPI.update(true);
            m_viewer->clearSelection();
            m_usedShapeConstraint = false;
        });
        /* end of modification */
        m_pdAPI.initConstraintsGUI();

        m_viewer->performLayout();
    }


    void updateGroupAvgPosition(){
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
                sum += m_pdAPI.getPositions().row(m_sortedIndex[j]);
            }
            m_avgPosition[i] = sum / num;
        }
    }


    /*begin of modifying*/
    void addEdgeSpringConstraintsSelection(ProjDyn::Scalar weight) {
        Surface_mesh* smesh = m_viewer->getMesh();
        std::vector<ProjDyn::ConstraintPtr> spring_constraints;
        const auto& selVerts = m_viewer->getSelectedVertices();

        for (auto edge : smesh->edges()) {
            // only add Edge Constraints to those edges whose two vertices are selected
            if (std::find(selVerts.begin(), selVerts.end(), smesh->vertex(edge, 0).idx()) == selVerts.end()) continue;
            if (std::find(selVerts.begin(), selVerts.end(), smesh->vertex(edge, 1).idx()) == selVerts.end()) continue;
            // The weight is set to the edge length
            ProjDyn::Scalar w = (m_initialPositions.row(smesh->vertex(edge, 0).idx()) - m_initialPositions.row(smesh->vertex(edge, 1).idx())).norm();
            if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                std::vector<Index> edge_inds;
                edge_inds.push_back(smesh->vertex(edge, 0).idx());
                edge_inds.push_back(smesh->vertex(edge, 1).idx());
                ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge_inds, w, m_initialPositions);
                spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
            }
        }
        m_pdAPI.addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Springs selection", spring_constraints, weight));
    }

    void addBaseAngleConstraintsSelection(ProjDyn::Scalar weight, bool right_side) {
        /***************************************************************************//**
        * A brief history of Doxygen-style banner comments.
        *
        * This is a Doxygen-style C-style "banner" comment. It starts with a "normal"
        * comment and is then converted to a "special" comment block near the end of
        * the first line. It is written this way to be more "visible" to developers
        * who are reading the source code.
        * This style of commenting behaves poorly with clang-format.
        *
        * @param theory Even if there is only one possible unified theory. it is just a
        *               set of rules and equations.
        ******************************************************************************/
        Surface_mesh* smesh = m_viewer->getMesh();
        static std::vector<ProjDyn::ConstraintPtr> angle_constraints;
        // The weight is set to the edge length
        ProjDyn::Scalar w = 1;
        std::vector<ProjDyn::Scalar> theta_target;
        theta_target.push_back(0.0);    // never used
        theta_target.push_back(0.0);
        theta_target.push_back(0.0);
        theta_target.push_back(0.0);
        theta_target.push_back(0.0);

        if(m_sortedIndex.size() == 0){
            // Default mode: if no vertices is selected by hand, then loading pre-selected vertices index from disk
            std::ifstream inFile;
            int x;
            inFile.open("../data/vertices_index.txt");
            if(!inFile){
                std::cout << "could not open/find file at ../data/vertices_index.txt" << std::endl;
                return;
            }
            for(int leg = 0; leg < 8; leg++){
                for(int i = 0; i < 4; i++){
                    inFile >> x;
                    m_sortedIndex.push_back(x);
                    m_groupNum.push_back(i+1);  //  first joint vertex is stored at 1st position of m_sortedIndex
                                                //  second joint is at 2nd position and so on
                }
                bool right = (leg > 3);         // right and left side of angle constrains will build opposite local coordinate
                ProjDyn::BaseAngleConstraint* esc = new ProjDyn::BaseAngleConstraint(m_sortedIndex, m_groupNum, w, m_initialPositions, right, 
                            theta_target, THETA_SCOPE);
                angle_constraints.push_back(std::shared_ptr<ProjDyn::BaseAngleConstraint>(esc));
                m_sortedIndex.clear();
                m_groupNum.clear();
                std::cout << leg + 1 << " legs has been added " << std::endl;
            }
        }
        else{
            // Custome mode: supporting select vertices of joint by hand
            ProjDyn::BaseAngleConstraint* esc = new ProjDyn::BaseAngleConstraint(m_sortedIndex, m_groupNum, w, m_initialPositions, right_side, 
                theta_target, THETA_SCOPE);
            angle_constraints.push_back(std::shared_ptr<ProjDyn::BaseAngleConstraint>(esc));
            std::cout << angle_constraints.size() << " legs has been added " << std::endl;
        }

        if(angle_constraints.size() == 8){
            m_pdAPI.addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Angle selection", angle_constraints, weight));
            std::cout << " Done !" << std::endl;
            angle_constraints.clear();
        }
    }
    /*end of modifying*/


    // Here, a flatness or smoothness constraint is added to the vertices
    // that have been selected in the viewer.
    // You can examine this method to understand how constraints can be created
    // and added to the viewer.
    void addFlatnessConstraintsToSelection(ProjDyn::Scalar weight) {
        if (!m_viewer) return;

        // List of indices of the vertices currently selected in the viewer
        const auto& selVerts = m_viewer->getSelectedVertices();

        // Vertex stars are simply a type of neighborhood information for each vertex.
        // Alternatively, Surface_mesh methods could be used here.
        // These have the advantage that we can pass each individual vertex' neighborhood
        // information to the constraint at a vertex, without passing the full mesh.
        auto vertexStars = makeVertexStars(m_pdAPI.getPositions().rows(), m_pdAPI.getTriangles());
        ProjDyn::Vector voronoiAreas = ProjDyn::vertexMasses(m_pdAPI.getPositions(), m_pdAPI.getTriangles());

        // If you want to work on the Surface_mesh to construct the constraint, it can be
        // received from the viewer:
        Surface_mesh* surf_mesh = m_viewer->getMesh();

        // Eigen matrices for positions and triangles of the mesh.
        // NOTE: these are the current, possibly deformed vertex positions from the simulator.
        // The original vertex positions are available as m_initialPositions (or from the
        // viewers Surface_mesh, but those need to be transformed to the simulator's format).
        auto pos = m_pdAPI.getPositions();
        auto tris = m_pdAPI.getTriangles();

        // Create a list of new constraints, one for each selected vertex/triangle
        std::vector<ConstraintPtr> newConstraints;
        for (auto v : surf_mesh->vertices()) {
            if (std::find(selVerts.begin(), selVerts.end(), v.idx()) == selVerts.end()) continue;
            if (surf_mesh->is_boundary(v)) continue;
            std::shared_ptr<FlatteningConstraint> flat_con = std::make_shared<FlatteningConstraint>
            (vertexStars[v.idx()], voronoiAreas(v.idx()) * 0.01, voronoiAreas(v.idx()), pos, tris);
            newConstraints.push_back(flat_con);
        }

        // Create a group from the list of constraints (which gets a name and a weight multiplier)
        auto con_group = std::make_shared<ConstraintGroup>("Flat", newConstraints, weight);

        // We use the ProjDynAPI to add this constraint to the simulator
        m_pdAPI.addConstraints(con_group);
    }

    // Add bending constraints to some vertices of the simulation
    void addIsometricOneRingConstraints( ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> iso_constraints;
        const ProjDyn::Positions& sim_verts = m_initialPositions;
        std::vector<ProjDyn::VertexStar> vStars = makeVertexStars(m_pdAPI.getPositions().rows(), m_pdAPI.getTriangles());
        ProjDyn::Vector voronoiAreas = ProjDyn::vertexMasses(m_initialPositions, m_pdAPI.getTriangles());
        Surface_mesh* smesh = m_viewer->getMesh();
        for (auto v : smesh->vertices()) {
            Index i = v.idx();
            if (smesh->is_boundary(v)) continue;
            if (i >= m_pdAPI.getNumVertices()) continue;
            // The weight is the voronoi area
            ProjDyn::Scalar w = voronoiAreas(i);
            if (w > 1e-6) {
                // The constraint is constructed, made into a shared pointer and appended to the list
                // of constraints.
                ProjDyn::IsometricOneRing* iso1 = new ProjDyn::IsometricOneRing(getStarIndices(vStars[i]), w, sim_verts);
                iso_constraints.push_back(std::shared_ptr<ProjDyn::IsometricOneRing>(iso1));
            }
        }
        m_pdAPI.addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Isometric 1-Rings", iso_constraints, weight));
    }

    // Turns the current shape into the rest shape and resets the constraints
    void makeCurrentShapeToRestShape() {
        Surface_mesh* viewer_mesh = m_viewer->getMesh();
        const Positions& cur_pos = m_pdAPI.getPositions();
        for (const auto& v : viewer_mesh->vertices()) {
            Index ind = v.idx();
            viewer_mesh->position(v).x = cur_pos(ind, 0);
            viewer_mesh->position(v).y = cur_pos(ind, 1);
            viewer_mesh->position(v).z = cur_pos(ind, 2);
        }
        setMesh();
    }

    // Adds a single position constraint that keeps the mesh centered.
    // The ShapeUpAPI is set up such that this constraint is always present if
    // there are no other position constraints, and disappears once there are.
    // Without this constraint, there is always a translational degree of freedom
    // left in the system, such that in practice, your mesh position will randomly
    // change during the local-global algorithm (the linear system is then indefinite).
    void addCenterConstraint() {
        const ProjDyn::Positions& pos = m_pdAPI.getPositions();
        ConstraintGroupPtr centerConGroup = std::make_shared<ConstraintGroup>
        ("Center", std::vector<ConstraintPtr>
            ({ std::make_shared<PositionConstraintGroup>
                (std::vector<Index>({ 0 }), 1., pos.row(0)) 
            }), 1.);
        m_pdAPI.addConstraints(centerConGroup);
    }

    // This method checks the "status" of each vertex for the visualization
    // in the viewer.
    // For now it simply marks vertices which take part in position constraints.
    void updateVertexStatus() {
        if (!m_viewer) return;

        Index numVerts = m_pdAPI.getPositions().rows();

        Eigen::Matrix<int, 1, -1> vStatus;
        vStatus.resize(1, numVerts);
        vStatus.setConstant(VertexStatus::Default);

        for (auto move_group : m_movableGroups) {
            for (auto move_con : move_group->constraints) {
                for (Index v : move_con->getIndices()) {
                    vStatus(0, v) = VertexStatus::Movable;
                }
            }
        }
        m_viewer->updateVertexStatus(vStatus);
    }

    // Create a position constraint group
    // This does not get added to the simulator right away!
    std::shared_ptr<PositionConstraintGroup> createPositionConstraint(const std::vector<Index>& indices, ProjDyn::Scalar weight) {
        const Positions& curPos = m_pdAPI.getPositions();
        Positions con_pos;
        con_pos.setZero(indices.size(), 3);
        Index ind = 0;
        for (Index v : indices) {
            con_pos.row(ind) = curPos.row(v);
            ind++;
        }
        // Create constraint groups and add them to the simulation
        std::shared_ptr<PositionConstraintGroup> con = std::shared_ptr<PositionConstraintGroup>(new PositionConstraintGroup(indices, weight, con_pos));
        return con;

    }

    // This creates and adds a position constraint group to the simulator.
    void addPositionConstraintGroup(const std::vector<Index>& vertInds) {
        std::shared_ptr<PositionConstraintGroup> con = createPositionConstraint(vertInds, SHAPEUP_MOVABLE_WEIGHT);
        auto conGroup = std::make_shared<ConstraintGroup>("Fixed Pos.", std::vector<ConstraintPtr>({ con }), 1.);
        m_pdAPI.addConstraints(conGroup);
        m_movableGroups.push_back(conGroup);
        // Remove center constraint, which is no longer required
        m_pdAPI.removeConstraints("Center");
        updateVertexStatus();        
    }

    // This method is called from the viewer, when the user alt + left-clicks on
    // some mesh vertices and drags them.
    // The 3d position of the mouse (this is defined in some clever way) is passed
    // along with the indices of the clicked vertices.
    // We check if the clicked vertices take part in any position constraint groups
    // and if so, change translate the target positions of the group with the most
    // clicked vertices.
    // The first time this event is called, a group is selected, and we keep translating
    // its target positions until the left mouse button is released (which calls the
    // method below).
    void mouseGrabEvent(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
        if (m_isUpdating) return;
        if (m_movableGroups.empty()) return;

        // Initialize grab if no group has been selected yet
        if (!m_isGrabbing) {
            // Compute the center of the grabbed vertices
            m_grabCenter.setZero();
            for (const Vector3f& gp : grabPos) {
                m_grabCenter(0) += gp(0);
                m_grabCenter(1) += gp(1);
                m_grabCenter(2) += gp(2);
            }
            m_grabCenter /= (ProjDyn::Scalar)grabPos.size();

            // Check which group of vertices has most grabbed vertices
            Index countCurrent = 0, countHighest = 0;
            Index indexHighest = 0, ind = 0;
            for (const auto& grp : m_movableGroups) {
                
                for (const auto& cur_con : grp->constraints) {
                    for (Index ind : grabbedVerts) {
                        if (std::find(cur_con->getIndices().begin(), cur_con->getIndices().end(), ind) != cur_con->getIndices().end()) {
                            countCurrent++;
                        }
                    }
                    if (countCurrent > countHighest) {
                        countHighest = countCurrent;
                        indexHighest = ind;
                    }
                    ind++;
                }
            }
            if (countHighest == 0) return;

            m_grabGroup = m_movableGroups[indexHighest];

            m_oldConstraintPos = std::dynamic_pointer_cast<PositionConstraintGroup>(m_grabGroup->constraints[0])->getTargetPositions();
            m_isGrabbing = true;
        }

        // If a group has been selected, translate its target positions
        if (m_grabGroup && m_isGrabbing) {
            m_isUpdating = true;
            // Update grab (translate the selected group by the amount that the center moved)
            Vector3 curGrabCenter(0., 0., 0.);
            for (const Vector3f& gp : grabPos) {
                curGrabCenter(0) += gp(0);
                curGrabCenter(1) += gp(1);
                curGrabCenter(2) += gp(2);
            }
            curGrabCenter /= (ProjDyn::Scalar)grabPos.size();
            Vector3 translate = curGrabCenter - m_grabCenter;
            Positions newTargetPos = m_oldConstraintPos;
            for (Index i = 0; i < newTargetPos.rows(); i++) {
                newTargetPos.row(i) += translate.transpose();
            }
            std::dynamic_pointer_cast<PositionConstraintGroup>(m_grabGroup->constraints[0])->setTargetPositions(newTargetPos);
            if (m_viewer) m_pdAPI.update(true);
            m_isUpdating = false;
        }
    }

    // Called when an alt + left-click event is finished. (See above.)
    void mouseGrabRelease() {
        m_isGrabbing = false;
        m_grabGroup = nullptr;
    }

    // When the simulator changes the positions using the local-global algorithm,
    // this is done internally, by modifying local n x 3 matrices containing the
    // positions.
    // To let the viewer know that those positions have been changed, we call the
    // method below.
    // This method does not actually change the viewer's Surface_mesh, but instead
    // directly uploads the updated positions to the GPU by passing them into an
    // OpenGL buffer from which the vertex positions are read.
    bool uploadPositions() {
        return m_pdAPI.uploadPositions();
    }

private:
    std::vector<ConstraintGroupPtr> m_movableGroups;

    bool m_isGrabbing = false;
    bool m_isUpdating = false;
    Index m_usedShapeConstraint = 0;
    ProjDyn::Vector3 m_grabCenter;
    Positions m_oldConstraintPos;
    Positions m_initialPositions;
    ConstraintGroupPtr m_grabGroup = nullptr;
    Viewer* m_viewer = nullptr;
    ProjDynAPI m_pdAPI;
    Button* m_buttonEdgeSprings = nullptr;
    bool m_isInitialized = false;
    /* begaining of self added varibles to record state */
    std::vector<Index> m_sortedIndex;               // store the vetices selected in sequence
    std::vector<int> m_groupNum;                    // store the counts of selected vertices at each joint
    std::vector<double> m_original_edge_length;     // store the original length of each limb part
    std::vector<ProjDyn::Vector3> m_avgPosition;    // store the average position of selected vertices at each joint
    /*end*/
};
