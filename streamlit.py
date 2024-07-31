import streamlit as st
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import pvlib
from pvlib.location import Location

# Streamlit UI for user inputs
st.title("Solar Panel Orientation Visualization")

latitude = st.number_input("Enter Latitude", value=11.279943232698631)
longitude = st.number_input("Enter Longitude", value=76.01286658771784)
timezone = 'Asia/Kolkata'  # Indian Standard Time (IST)
start_date = st.date_input("Start Date", value=pd.to_datetime("2024-07-27"))
end_date = st.date_input("End Date", value=pd.to_datetime("2024-07-29"))

# Add a button to trigger the calculation
if st.button("Calculate"):
    # Convert dates to strings for use in calculations
    start_date = start_date.strftime('%Y-%m-%d')
    end_date = end_date.strftime('%Y-%m-%d')

    # Create a Location object
    location = Location(latitude, longitude, timezone)

    # Generate the solar position data for the specified date range
    times = pd.date_range(start=start_date, end=end_date, freq='1min', tz=timezone)
    solar_position = location.get_solarposition(times)

    # Calculate the daily optimal azimuth angles
    daily_optimal_azimuths = solar_position.groupby(solar_position.index.date)['elevation'].idxmax().apply(
        lambda x: solar_position.loc[x, 'azimuth'])

    # Calculate the average optimal azimuth angle for the period
    average_optimal_azimuth = daily_optimal_azimuths.mean()

    # Calculate the compromise azimuth
    solar_noon_times = solar_position.between_time('11:00', '13:00')
    compromise_azimuth = solar_noon_times['azimuth'].mean()

    # Set the panel direction to this compromise azimuth angle
    panel_direction = compromise_azimuth

    # Rotation matrix to align the prism with the panel direction
    def rotation_matrix(angle_degrees):
        angle_radians = np.radians(angle_degrees)
        cos_angle = np.cos(angle_radians)
        sin_angle = np.sin(angle_radians)
        return np.array([
            [cos_angle, -sin_angle, 0],
            [sin_angle, cos_angle, 0],
            [0, 0, 1]
        ])

    # Apply rotation to vertices
    def rotate_vertices(vertices, angle_degrees):
        R = rotation_matrix(angle_degrees)
        return np.dot(vertices, R.T)

    # Define the triangular prism vertices and faces
    base = 0.065  # Reduced by 10 times
    hypotenuse = 0.365  # Reduced by 10 times
    angle = 79  # base angle in degrees
    height = np.tan(np.radians(angle)) * base

    # Vertices of the triangular prism centered around the origin
    vertices = np.array([
        [0, -0.05, 0],  # A
        [base, -0.05, 0],  # B
        [0, -0.05, height],  # C
        [0, 0.05, 0],  # D
        [base, 0.05, 0],  # E
        [0, 0.05, height],  # F
    ])

    # Rotate the vertices to align with the panel direction
    vertices = rotate_vertices(vertices, -panel_direction + 90)  # Adjusting rotation

    x = []
    y = []
    z = []
    i = []
    j = []
    k = []

    def add_face(face_vertices):
        face_vertices = np.array(face_vertices)
        x.extend(face_vertices[:, 0])
        y.extend(face_vertices[:, 1])
        z.extend(face_vertices[:, 2])
        offset = len(x) - 4
        i.extend([offset, offset + 1, offset + 2, offset + 3])
        j.extend([offset + 1, offset + 2, offset + 3, offset])
        k.extend([offset + 2, offset + 3, offset, offset + 1])

    # Define the rectangular faces of the prism
    faces = [
        [vertices[0], vertices[1], vertices[4], vertices[3]],  # ABED
        [vertices[1], vertices[2], vertices[5], vertices[4]],  # BCEF
        [vertices[2], vertices[0], vertices[3], vertices[5]],  # CADF
    ]

    # Add the rectangular faces to the plot
    for face in faces:
        add_face(face)

    # Also add triangular faces
    x_tri = []
    y_tri = []
    z_tri = []
    i_tri = []
    j_tri = []
    k_tri = []

    def add_triangular_face(face_vertices):
        face_vertices = np.array(face_vertices)
        x_tri.extend(face_vertices[:, 0])
        y_tri.extend(face_vertices[:, 1])
        z_tri.extend(face_vertices[:, 2])
        offset = len(x_tri) - 3
        i_tri.extend([offset, offset + 1, offset + 2])
        j_tri.extend([offset + 1, offset + 2, offset])
        k_tri.extend([offset + 2, offset, offset + 1])

    add_triangular_face([vertices[0], vertices[1], vertices[2]])  # ABC
    add_triangular_face([vertices[3], vertices[4], vertices[5]])  # DEF

    fig_3d = go.Figure()

    fig_3d.add_trace(go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        facecolor=['lightblue'] * len(i),
        opacity=1.0,
        name='Triangular Prism Rectangles'
    ))

    fig_3d.add_trace(go.Mesh3d(
        x=x_tri,
        y=y_tri,
        z=z_tri,
        i=i_tri,
        j=j_tri,
        k=k_tri,
        facecolor=['lightblue'] * len(i_tri),
        opacity=1.0,
        name='Triangular Prism Triangles'
    ))

    # Draw the device marker at the origin
    fig_3d.add_trace(go.Scatter3d(
        x=[0],
        y=[0],
        z=[0],
        mode='markers+text',
        text=['Device\n(center)'],
        textposition='top center',
        marker=dict(size=12, color='blue'),
        name='Device Center'
    ))

    # Add cardinal directions
    cardinals = {'N': (0, 1), 'S': (0, -1), 'E': (1, 0), 'W': (-1, 0)}
    for direction, (x, y) in cardinals.items():
        fig_3d.add_trace(go.Scatter3d(
            x=[x],
            y=[y],
            z=[0],
            mode='text',
            text=[direction],
            textposition='middle center',
            showlegend=False
        ))

    # Correct the sun path plotting
    sun_path_times = pd.date_range(start=start_date, end=end_date, freq='H', tz=timezone)
    sun_positions = location.get_solarposition(sun_path_times)

    # Group by date and plot the sun path for each day
    for date, group in sun_positions.groupby(sun_positions.index.date):
        fig_3d.add_trace(go.Scatter3d(
            x=0.5 * np.sin(np.radians(group['azimuth'])),  # sin for x (east-west)
            y=0.5 * np.cos(np.radians(group['azimuth'])),  # cos for y (north-south)
            z=group['elevation'] / 90.0,  # Normalize elevation to fit in the plot
            mode='lines',
            line=dict(width=2),
            name=f'Sun Path on {date}'
        ))

    # Set 3D plot layout without axes
    fig_3d.update_layout(
        title=f"Optimal Solar Panel Orientation<br>Panel Direction: {panel_direction:.2f} degrees from North",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectratio=dict(x=1, y=1, z=0.5),
        ),
        showlegend=True,
        autosize=False,
        width=800,
        height=800
    )

    # Display the 3D plot in Streamlit
    st.plotly_chart(fig_3d)

    # Create 2D top view plot
    fig_2d = go.Figure()

    # Draw the device marker at the origin in 2D
    fig_2d.add_trace(go.Scatter(
        x=[0],
        y=[0],
        mode='markers+text',
        text=['Device\n(center)'],
        textposition='top center',
        marker=dict(size=12, color='blue'),
        name='Device Center'
    ))

    # Add cardinal directions in 2D
    for direction, (x, y) in cardinals.items():
        fig_2d.add_trace(go.Scatter(
            x=[x],
            y=[y],
            mode='text+markers',
            text=[direction],
            textposition='middle center',
            showlegend=False,
            marker=dict(size=8, color='black', symbol='arrow-bar-up')
        ))

    # Add North-South and East-West lines
    fig_2d.add_shape(type="line", x0=0, y0=1, x1=0, y1=-1, line=dict(color="black", width=1))
    fig_2d.add_shape(type="line", x0=-1, y0=0, x1=1, y1=0, line=dict(color="black", width=1))

    # Add the rectangular faces of the prism to the 2D plot
    for face in faces:
        x = [vertex[0] for vertex in face]
        y = [vertex[1] for vertex in face]
        x.append(face[0][0])
        y.append(face[0][1])
        fig_2d.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines',
            fill='toself',
            fillcolor='lightblue',
            line=dict(color='black'),
            name='Triangular Prism Face'
        ))

    # Mark the direction and angle in the 2D plot
    fig_2d.add_shape(type="line",
                     x0=0, y0=0, x1=np.sin(np.radians(panel_direction)), y1=np.cos(np.radians(panel_direction)),
                     line=dict(color="red", width=3))

    fig_2d.add_annotation(
        x=np.sin(np.radians(panel_direction)),
        y=np.cos(np.radians(panel_direction)),
        text=f"{panel_direction:.2f}°",
        showarrow=True,
        arrowhead=2
    )

    # Set 2D plot layout
    fig_2d.update_layout(
        title=f"2D Top View of Triangular Prism<br>Panel Direction: {panel_direction:.2f} degrees from North",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        autosize=False,
        width=800,
        height=800
    )

    # Display the 2D plot in Streamlit
    st.plotly_chart(fig_2d)
