# digitalGeometryProcessing
- The project of this course is about making an extension of lab's delivered computer graphics framework Projective Dynamics and ShapeOp.
- My proposal is to implement a realistic kinematics limb movement simulation on a spider mesh.
- In addition, I realized a simple walking gate controlled by cosine signals, the computing efficiency is high enough
to fufill the need of real-time rendering.

In short, the magic is to animate a static character, only with input of its 3D mesh. But pitifully, my current implementaion hardcode the algorithm with input mesh data. So it only works for my spider mesh. It's interesting and challenging to generalize its ability to any input mesh, with minimum manual setup. Maybe need some machine learning technique, up to your imagination. 

Let's look at some demos of my project!

!(spider1)[https://firebasestorage.googleapis.com/v0/b/steam-key-269816.appspot.com/o/Enregistrement%20de%20l%E2%80%99e%CC%81cran%202019-12-16%20a%CC%80%2013.20.09.gif?alt=media&token=2c71c105-5550-48ee-b580-4a60b079b3a8]
!(spider2)[https://firebasestorage.googleapis.com/v0/b/steam-key-269816.appspot.com/o/Enregistrement%20de%20l%E2%80%99e%CC%81cran%202019-12-16%20a%CC%80%2013.07.15.gif?alt=media&token=d30e6c68-0928-48c9-b0fb-cbbea68ed15d]
